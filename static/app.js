const dz = document.getElementById('dropzone');
const form = document.getElementById('job-form');

// helper: bind change handler to the current file input element
function bindFileInput() {
  const input = document.getElementById('file-input');
  if (!input) return;
  // avoid double-binding
  if (input._bound) return;
  input.addEventListener('change', e => {
    if (input.files && input.files[0]) {
      handleFile(input.files[0]);
    }
  });
  input._bound = true;
}

let isUploading = false;
let stagedRunId = null;
let stagedFilename = null;

function setDropzonePreview(file) {
  if (!dz) return;
  dz.classList.add('has-file');
  dz.innerHTML = `
    <div class="file-info">
      <strong>${file.name}</strong>
      <div class="muted">${(file.size / (1024*1024)).toFixed(2)} MB</div>
      <div class="muted">Click to replace</div>
    </div>
    <input id="file-input" type="file" name="vcf" accept=".vcf,.gz,.bcf" style="display:none" />
  `;
  // ensure the new input is bound
  bindFileInput();
}

function resetDropzone() {
  if (!dz) return;
  dz.classList.remove('has-file');
  dz.innerHTML = `
    <p>Drag & drop your <strong>.vcf</strong>, <strong>.vcf.gz</strong> or <strong>.bcf</strong> here</p>
    <p class="muted">or click to browse…</p>
    <input id="file-input" type="file" name="vcf" accept=".vcf,.gz,.bcf" required />
  `;
  // Bind change handler to the newly created input
  bindFileInput();
}

function updateUploadStatus(file, loaded, total) {
  const status = document.getElementById('upload-status');
  const bar = document.getElementById('upload-bar');
  const pct = document.getElementById('upload-percent');
  const fname = document.getElementById('upload-filename');
  const bytes = document.getElementById('upload-bytes');
  const submitBtn = document.getElementById('submit-btn');

  if (status) status.hidden = false;
  if (fname) fname.textContent = file.name;
  if (pct && typeof loaded === 'number' && typeof total === 'number') {
    const percent = Math.round((loaded / total) * 100);
    pct.textContent = percent + '%';
    if (bar) bar.style.width = percent + '%';
    if (bytes) bytes.textContent = (loaded / (1024*1024)).toFixed(1) + ' / ' + (total / (1024*1024)).toFixed(1) + ' MB';
  }
  if (submitBtn) submitBtn.disabled = isUploading;
}

function handleFile(file) {
  if (!file) return;
  // show preview and create a fresh hidden input
  setDropzonePreview(file);
  // assign the file to the newly-created input element so it's available for FormData
  try {
    const cur = document.getElementById('file-input');
    if (cur) {
      const dataTransfer = new DataTransfer();
      dataTransfer.items.add(file);
      cur.files = dataTransfer.files;
    }
  } catch (e) {
    // older browsers: nothing to do
  }

  // Immediately upload after validating required selects (species & chromosome)
  const species = document.getElementById('species-select') ? document.getElementById('species-select').value : (form.querySelector('select[name=species]') && form.querySelector('select[name=species]').value);
  const chromosome = document.getElementById('chromosome-select') ? document.getElementById('chromosome-select').value : (form.querySelector('select[name=chromosome]') && form.querySelector('select[name=chromosome]').value);

  // Stage the file immediately (send to /upload) using XHR so we can show progress.
  const fd = new FormData();
  fd.append('vcf', file);
  isUploading = true;
  updateUploadStatus(file, 0, file.size);

  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/upload', true);
  xhr.responseType = 'json';
  xhr.upload.onprogress = function (ev) {
    if (ev.lengthComputable) {
      updateUploadStatus(file, ev.loaded, ev.total);
    }
  };
  xhr.onerror = function () {
    isUploading = false;
    alert('Upload failed. Please try again.');
    resetDropzone();
  };
  xhr.onload = function () {
    isUploading = false;
    try {
      const data = xhr.response;
      if (!data && xhr.responseText) {
        // some browsers don't parse json when responseType=json and status != 200
        try { data = JSON.parse(xhr.responseText); } catch (e) { data = null; }
      }
      if (data && data.run_id && data.filename) {
        stagedRunId = data.run_id;
        stagedFilename = data.filename;
        const fnameEl = document.getElementById('upload-filename');
        if (fnameEl) fnameEl.textContent = stagedFilename + ' (staged)';
        updateUploadStatus(file, file.size, file.size);
      } else if (data && data.error) {
        alert('Upload failed: ' + data.error);
        resetDropzone();
      } else if (xhr.status < 200 || xhr.status >= 300) {
        alert('Upload failed: server returned ' + xhr.status);
        resetDropzone();
      }
    } catch (e) {
      alert('Upload failed: ' + e);
      resetDropzone();
    }
  };
  xhr.send(fd);
}

function uploadFile(file) {
  if (!form || !file) return;
  if (isUploading) return; // avoid duplicate uploads
  isUploading = true;

  const species = form.querySelector('select[name=species]').value;
  const chromosome = form.querySelector('select[name=chromosome]').value;
  const grid = form.querySelector('input[name=grid]').value || '300';

  updateUploadStatus(file, 0, file.size);

  const data = new FormData();
  data.append('species', species);
  data.append('chromosome', chromosome);
  data.append('grid', grid);
  data.append('vcf', file);

  const xhr = new XMLHttpRequest();
  xhr.open('POST', form.action, true);
  xhr.responseType = 'document';
  xhr.upload.onprogress = function (ev) {
    if (ev.lengthComputable) {
      updateUploadStatus(file, ev.loaded, ev.total);
    }
  };
  xhr.onerror = function () {
    isUploading = false;
    updateUploadStatus(file, 0, file.size);
    alert('Upload failed. Please try again.');
  };
  xhr.onload = function () {
    isUploading = false;
    updateUploadStatus(file, file.size, file.size);
    if (xhr.status >= 200 && xhr.status < 300) {
      // Replace the page with the server-rendered result
      document.open();
      document.write(xhr.response.documentElement.outerHTML);
      document.close();
    } else {
      alert('Server error: ' + xhr.status);
    }
  };
  xhr.send(data);
}

if (dz) {
  // drag visual state
  ['dragenter','dragover'].forEach(evt => dz.addEventListener(evt, e => { e.preventDefault(); e.stopPropagation(); dz.classList.add('drag'); }));
  ['dragleave','drop'].forEach(evt => dz.addEventListener(evt, e => { e.preventDefault(); e.stopPropagation(); dz.classList.remove('drag'); }));
  dz.addEventListener('drop', e => {
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      // assign files and immediately handle
      handleFile(e.dataTransfer.files[0]);
    }
  });

  // if user clicks the dropzone to open file picker, keep default input behavior
  dz.addEventListener('click', (e) => {
    // If the user clicked the actual file input (or an element inside it),
    // don't re-open the picker — that causes the dialog to appear twice.
    if (e && e.target) {
      try {
        if (e.target.closest && e.target.closest('#file-input')) {
          return;
        }
      } catch (err) {
        // ignore and continue
      }
    }
    const pick = document.getElementById('file-input');
    if (pick) pick.click();
  });

  // bind change on the input (ensures we bind the current element)
  bindFileInput();

  // Fetch chromosomes when species is changed
  const speciesSelect = document.getElementById('species-select') || document.querySelector('select[name=species]');
  const chromSelect = document.getElementById('chromosome-select');
  if (speciesSelect && chromSelect) {
    function updateChromosomes() {
      const species = speciesSelect.value;
      // don't fetch when placeholder is selected
      if (!species) {
        chromSelect.innerHTML = '<option value="">Select a chromosome</option>';
        return;
      }
      fetch(`/chromosomes?species=${encodeURIComponent(species)}`)
        .then(resp => resp.json())
        .then(data => {
          if (data.chromosomes && Array.isArray(data.chromosomes)) {
            chromSelect.innerHTML = '<option value="">Select a chromosome</option>';
            data.chromosomes.forEach(chrom => {
              chromSelect.innerHTML += `<option value="${chrom}">${chrom}</option>`;
            });
          } else {
            chromSelect.innerHTML = '<option value="">No chromosomes found</option>';
          }
        })
        .catch(() => {
          chromSelect.innerHTML = '<option value="">Error loading chromosomes</option>';
        });
    }
    speciesSelect.addEventListener('change', updateChromosomes);
    // Initial population on page load
    updateChromosomes();
  }
}

// Keep submit behavior: if user presses Analyze and a file is present, start upload (or ignore if already uploading)
if (form) {
  form.addEventListener('submit', function (e) {
    const input = document.getElementById('file-input');
    const species = document.getElementById('species-select') ? document.getElementById('species-select').value : (form.querySelector('select[name=species]') && form.querySelector('select[name=species]').value);
    const chromosome = document.getElementById('chromosome-select') ? document.getElementById('chromosome-select').value : (form.querySelector('select[name=chromosome]') && form.querySelector('select[name=chromosome]').value);

    const missing = [];
    if (!species) missing.push('species');
    if (!chromosome) missing.push('chromosome');

    if (missing.length > 0) {
      e.preventDefault();
      if (missing.length === 2) {
        alert('Please select a species and a chromosome before analyzing.');
      } else if (missing[0] === 'species') {
        alert('Please select a species before analyzing.');
      } else {
        alert('Please select a chromosome before analyzing.');
      }
      return;
    }

    // If a staged upload exists, submit the analyze form including run_id & filename
    if (stagedRunId && stagedFilename) {
      e.preventDefault();
      // build a FormData to post to /analyze using staged file
      const fd = new FormData();
      fd.append('run_id', stagedRunId);
      fd.append('uploaded_filename', stagedFilename);
      fd.append('species', species);
      fd.append('chromosome', chromosome);
      const grid = form.querySelector('input[name=grid]').value || '300';
      fd.append('grid', grid);

      // Start polling the run log while analysis runs, then POST to /analyze.
      e.preventDefault();
      // ensure analysis-status container exists
      let status = document.getElementById('analysis-status');
      if (!status) {
        status = document.createElement('div');
        status.id = 'analysis-status';
        status.style = 'white-space:pre-wrap; background:#111; color:#efe; padding:10px; margin-top:12px; max-height:300px; overflow:auto; font-family: monospace;';
        form.parentNode.insertBefore(status, form.nextSibling);
      }
      status.textContent = 'Starting analysis...\n';

      let poller = null;
      function startPolling(runId) {
        if (poller) return;
        poller = setInterval(() => {
          fetch(`/runs/${encodeURIComponent(runId)}/tail`)
            .then(r => r.json())
            .then(data => {
              if (data && data.lines !== undefined) {
                status.textContent = data.lines;
                // scroll to bottom
                status.scrollTop = status.scrollHeight;
              }
              if (data && data.done) {
                // let the analyze response replace the page; stop polling
                clearInterval(poller);
                poller = null;
              }
            })
            .catch(() => {
              // ignore polling errors; keep trying
            });
        }, 1000);
      }

      // start polling using the staged run id
      if (stagedRunId) startPolling(stagedRunId);

      // POST to /analyze which now starts the background job and returns JSON
      fetch(form.action, { method: 'POST', body: fd })
        .then(r => {
          if (!r.ok) throw new Error('Server returned ' + r.status);
          return r.json();
        })
        .then(data => {
          if (!data || !data.run_id) throw new Error('Invalid response from server');
          const runId = data.run_id;

          // ensure polling is running
          if (stagedRunId) stagedRunId = runId; else stagedRunId = runId;

          // continue polling for logs until done, then fetch final HTML from /runs/<run_id>/final
          let retries = 0;
          const maxRetries = 5;
          function fetchFinal() {
            fetch(`/runs/${encodeURIComponent(runId)}/final`)
              .then(r => {
                if (r.status === 202) {
                  // not ready yet, wait and retry
                  return null;
                }
                if (!r.ok) throw new Error('Server returned ' + r.status);
                return r.text();
              })
              .then(html => {
                if (html) {
                  if (poller) { clearInterval(poller); poller = null; }
                  document.open();
                  document.write(html);
                  document.close();
                }
              })
              .catch(err => {
                // Network errors can happen; retry a few times before failing
                retries += 1;
                if (retries <= maxRetries) {
                  setTimeout(fetchFinal, 1000 * retries);
                } else {
                  if (poller) { clearInterval(poller); poller = null; }
                  alert('Analysis failed: ' + err + '. Check your network or server logs.');
                }
              });
          }

          // start a second-level poll that tries to fetch the final page every 2s
          const finalPoll = setInterval(() => fetchFinal(), 2000);

          // also ensure we stop the finalPoll when poller indicates done via /runs/<run_id>/tail
          const origPoller = poller;
          // fetchFinal will clear poller and replace the page when ready
        })
        .catch(err => {
          if (poller) { clearInterval(poller); poller = null; }
          alert('Analysis failed: ' + err);
        });
      return;
    }

    // if no staged file, but a file-selected input exists, let upload happen via existing uploadFile
    if (!input || !input.files || input.files.length === 0) return;
    e.preventDefault();
    const file = input.files[0];
    if (!isUploading) uploadFile(file);
  });
}