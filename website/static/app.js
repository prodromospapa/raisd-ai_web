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
// Track active upload XHR so we can abort when user replaces file
let lastUploadXhr = null;
// Track analysis polling interval so we can stop it when user replaces staged run
let analysisPoller = null;
// remember the last server-suggested chromosome so we can preserve it when
// the chromosome list is refreshed after species selection
let lastSuggestedChromosome = null;
// suppress showing the auto-selection popover while we refresh the
// chromosome list (e.g. after species change)
let suppressAutoMessageDuringRefresh = false;

// ensure the CSS for the chromosome glow animation is injected once
function ensureChromGlowCSS() {
  if (document.getElementById('chrom-glow-css')) return;
  const css = `
  /* Neon-style autofill pulse + halo inspired by the provided image.
     Uses a subtle scale + halo fade to draw attention without blocking UI. */
  @keyframes neonPulse {
    0% {
      opacity: 0; transform: scale(0.995);
      box-shadow: 0 0 0 rgba(0,180,255,0);
    }
    30% {
      opacity: 1; transform: scale(1.01);
      box-shadow: 0 0 24px 6px rgba(0,180,255,0.28);
    }
    60% {
      opacity: 0.9; transform: scale(1.006);
      box-shadow: 0 0 18px 4px rgba(0,140,220,0.22);
    }
    100% {
      opacity: 0; transform: scale(0.995);
      box-shadow: 0 0 0 rgba(0,180,255,0);
    }
  }

  .chrom-glow-anim { position: relative; display: inline-block; z-index: 0; }

  /* inner rim for crisp neon edge */
  .chrom-glow-anim::before {
    content: ""; position: absolute; inset: 0; pointer-events: none; border-radius: 10px;
    border: 1px solid rgba(0,200,255,0.35);
    mix-blend-mode: screen; opacity: 0; transform: scale(0.998);
  }

  /* outer halo that pulses and fades */
  .chrom-glow-anim::after {
    content: ""; position: absolute; left: -8px; right: -8px; top: -6px; bottom: -6px; pointer-events: none;
    border-radius: 14px; opacity: 0; transform: scale(0.996);
    background: radial-gradient(closest-side, rgba(0,200,255,0.16), rgba(0,100,180,0.03));
    filter: blur(6px);
  }

  /* apply the neon pulse on run */
  .chrom-glow-anim.chrom-glow-run::before { animation: neonPulse 900ms cubic-bezier(.2,1,.3,1) forwards; }
  .chrom-glow-anim.chrom-glow-run::after  { animation: neonPulse 1100ms cubic-bezier(.2,1,.3,1) forwards; }

  /* subtle text glow for controls with text when active */
  .chrom-glow-anim.chrom-glow-run select, .chrom-glow-anim.chrom-glow-run input { text-shadow: 0 0 8px rgba(0,180,255,0.15); }

  @media (prefers-reduced-motion: reduce) {
    .chrom-glow-anim::before, .chrom-glow-anim::after { animation: none !important; opacity: 0 !important; filter: none !important; }
  }
  `;
  const style = document.createElement('style');
  style.id = 'chrom-glow-css';
  style.appendChild(document.createTextNode(css));
  (document.head || document.documentElement).appendChild(style);
}

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

// Show a short, native validation-style message on a select element.
// This creates a temporary customValidity and calls reportValidity so browsers
// display the validation tooltip similarly to the screenshot.
// animate the chromosome select to show a brief green highlight when
// a suggested chromosome is applied. No popover or message is shown.
function animateChromosomeSelect(selectEl) {
  if (!selectEl) return;
  if (suppressAutoMessageDuringRefresh) return;
  try {
    ensureChromGlowCSS();
    // Prefer a dedicated wrapper if present, otherwise fall back to closest label/parent
    const wrapper = selectEl.parentElement && selectEl.parentElement.classList && selectEl.parentElement.classList.contains('chromosome-list-wrapper') ? selectEl.parentElement : null;
    const container = wrapper || (typeof selectEl.closest === 'function' && selectEl.closest('label')) || selectEl.parentElement || selectEl;
    if (!container) return;
    // ensure base class present on container, then trigger run class which animates ::after
    container.classList.add('chrom-glow-anim');

    // trigger glow-only animation on next frame
    requestAnimationFrame(() => {
      container.classList.add('chrom-glow-run');
      const cleanup = () => { try { container.classList.remove('chrom-glow-run'); } catch(_){} };
      const onAnimEnd = () => { cleanup(); container.removeEventListener('animationend', onAnimEnd); };
      container.addEventListener('animationend', onAnimEnd);
  // safety fallback (longer than the animation durations)
  const tid = setTimeout(() => { try { cleanup(); container.removeEventListener('animationend', onAnimEnd); } catch(_){}; clearTimeout(tid); }, 1400);
    });
  } catch (e) { /* ignore */ }
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
    if (bar) {
      // Keep animation active even at 100% until XHR onload/onerror explicitly removes it.
      if (!bar.classList.contains('animate')) {
        bar.classList.add('animate');
      }
    }
  }
  if (submitBtn) submitBtn.disabled = isUploading;
}

function handleFile(file) {
  if (!file) return;
  // If there's an in-flight upload, abort it — the new file replaces the old
  try {
    if (lastUploadXhr && typeof lastUploadXhr.abort === 'function') {
      try { lastUploadXhr.abort(); } catch(_){}
      lastUploadXhr = null;
    }
  } catch(_){}
  // If a previous staged run exists, ask server to clean it up (fire-and-forget)
  try {
    if (stagedRunId) {
      fetch(`/runs/${encodeURIComponent(stagedRunId)}/cleanup`, { method: 'POST', keepalive: true }).catch(()=>{});
      stagedRunId = null;
      stagedFilename = null;
      // Clear any analysis poller that may be watching the previous run
      try { if (analysisPoller) { clearInterval(analysisPoller); analysisPoller = null; } } catch(_){}
      // Reset UI upload status briefly
      try { const fnameEl = document.getElementById('upload-filename'); if (fnameEl) fnameEl.textContent = ''; } catch(_){}
    }
  } catch(_){}
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
  try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.add('animate'); } catch(_) {}

  const xhr = new XMLHttpRequest();
  lastUploadXhr = xhr;
  xhr.open('POST', '/upload', true);
  xhr.responseType = 'json';
  xhr.upload.onprogress = function (ev) {
    if (ev.lengthComputable) {
      updateUploadStatus(file, ev.loaded, ev.total);
    }
  };
  xhr.onerror = function () {
    isUploading = false;
    lastUploadXhr = null;
    alert('Upload failed. Please try again.');
    resetDropzone();
    try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
  };
  xhr.onload = function () {
  isUploading = false;
    lastUploadXhr = null;
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
        try {
          const bar = document.getElementById('upload-bar');
          if (bar) {
            bar.classList.remove('animate');
            bar.classList.add('completed');
            setTimeout(()=>{ bar.classList.remove('completed'); }, 1200);
          }
        } catch(_) {}
        // If server returned contigs, populate chromosome selector and suggest one
        try {
          if (data.contigs && Array.isArray(data.contigs) && data.contigs.length > 0) {
            const chromSelect = document.getElementById('chromosome-select');
            // Build a display set: prefer stripped form (drop leading 'chr') for the UI
            const displaySet = new Set();
            data.contigs.forEach(c0 => {
              const c = String(c0);
              // prefer the non-'chr' form for selection/display
              const display = /^chr/i.test(c) ? c.replace(/^chr/i, '') : c;
              displaySet.add(display);
            });
            const normalized = Array.from(displaySet);
            if (chromSelect) {
              // Put a detection message in the placeholder option so it appears in the select
              const shown = normalized.slice(0, 20).join(', ') + (normalized.length > 20 ? ', ...' : '');
              // normalize suggested chromosome for display as well
              const sugNorm = data.suggested_chromosome ? String(data.suggested_chromosome).replace(/^chr/i, '') : null;
              const placeholderText = 'Detected: ' + shown + (sugNorm ? (' — suggested: ' + String(sugNorm)) : '');
                // Use selected+hidden on the placeholder so it displays in the closed select
                // but doesn't appear in the open dropdown list.
                chromSelect.innerHTML = `<option value="" selected hidden>${placeholderText}</option>`;
              normalized.forEach(v => {
                const opt = document.createElement('option');
                opt.value = v;
                opt.textContent = v;
                chromSelect.appendChild(opt);
              });
              // If server suggested a chromosome, select it if present (use normalized form)
              if (sugNorm) {
                try { chromSelect.value = sugNorm; } catch(_) { }
                // animate the chromosome select to indicate the selection
                try { animateChromosomeSelect(chromSelect); } catch(_) {}
                // remember this suggestion so species->chromosome refresh won't clobber it
                try { lastSuggestedChromosome = String(sugNorm); } catch (_) {}
              }
            }

            // If the server suggested a chromosome, set it and animate the select so
            // the user can confirm before starting analysis. Do NOT auto-submit;
            // the user must click "Analyze" to begin processing.
            try {
              if (data.suggested_chromosome) {
                const chromSel = document.getElementById('chromosome-select');
                if (chromSel) {
                  chromSel.value = String(data.suggested_chromosome);
                  try { animateChromosomeSelect(chromSel); } catch(_) {}
                  try { lastSuggestedChromosome = String(data.suggested_chromosome); } catch(_) {}
                }
              }
            } catch (e) { /* ignore */ }
          }
        } catch (e) {
          // ignore UI update errors
          console.warn('Failed to apply contig suggestion:', e);
        }
      } else if (data && data.error) {
  alert('Upload failed: ' + data.error);
        resetDropzone();
        try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
      } else if (xhr.status < 200 || xhr.status >= 300) {
        alert('Upload failed: server returned ' + xhr.status);
        resetDropzone();
        try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
      }
    } catch (e) {
      alert('Upload failed: ' + e);
      resetDropzone();
      try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
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
  xhr.responseType = 'json';
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
    try {
      let data = xhr.response;
      if (!data && xhr.responseText) { try { data = JSON.parse(xhr.responseText); } catch(_){} }
      if (xhr.status >= 200 && xhr.status < 300 && data && (data.wait_url || data.run_id)) {
        const next = data.wait_url || (`/runs/${encodeURIComponent(data.run_id)}/final`);
        window.location.href = next;
      } else if (data && data.error) {
        alert('Server error: ' + data.error);
      } else {
        alert('Server error: ' + xhr.status);
      }
    } catch (e) {
      alert('Server error: ' + e);
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
    // Track whether the user manually selected a chromosome so we don't clobber it
    try {
      chromSelect._userSelected = false;
      chromSelect.addEventListener('change', () => { try { chromSelect._userSelected = true; lastSuggestedChromosome = null; } catch(_){} });
    } catch (e) { /* ignore */ }
    // Helper: populate chromSelect options array. If selectSuggestion is true and
    // a server suggestion exists (lastSuggestedChromosome), select it when present.
    function populateChromOptions(chromosomes, selectSuggestion) {
      try {
        // remember current selection so we can preserve it when rebuilding options
        const currentSelection = chromSelect ? String(chromSelect.value || '') : '';
  chromSelect.innerHTML = '<option value="" selected hidden>Select a chromosome</option>';
        const disp = Array.from(new Set(chromosomes.map(c => String(c).replace(/^chr/i, ''))));
        disp.forEach(chrom => {
          const opt = document.createElement('option');
          opt.value = chrom;
          opt.textContent = chrom;
          chromSelect.appendChild(opt);
        });
        // Preserve user's existing choice only if the user actually picked it.
        // If the selection was previously auto-suggested, prefer the new server suggestion.
        if (currentSelection && chromSelect._userSelected) {
          const foundOld = Array.from(chromSelect.options).some(o => o.value === currentSelection);
          if (foundOld) {
            chromSelect.value = currentSelection;
            return;
          }
        }
        // If no prior selection to preserve, optionally apply server suggestion
        if (selectSuggestion && lastSuggestedChromosome) {
          try {
            const normSug = String(lastSuggestedChromosome).replace(/^chr/i, '');
            const found = Array.from(chromSelect.options).some(o => o.value === normSug);
            if (found) chromSelect.value = normSug;
          } catch (e) { /* ignore */ }
        }
      } catch (e) {
        chromSelect.innerHTML = '<option value="" selected hidden>Error building chromosome list</option>';
      }
    }

    function updateChromosomes() {
      const species = speciesSelect.value;
      // don't fetch when placeholder is selected
      if (!species) {
        chromSelect.innerHTML = '<option value="" selected hidden>Select a chromosome</option>';
        return;
      }
      // prevent showing the auto-selection popover while we refresh options
      suppressAutoMessageDuringRefresh = true;
      fetch(`/chromosomes?species=${encodeURIComponent(species)}`)
        .then(resp => resp.json())
        .then(data => {
          // done fetching: allow messages again after we processed options
          suppressAutoMessageDuringRefresh = false;
          if (data.chromosomes && Array.isArray(data.chromosomes)) {
            populateChromOptions(data.chromosomes, true);
          } else {
            chromSelect.innerHTML = '<option value="" selected hidden>No chromosomes found</option>';
          }
        })
        .catch(() => {
          suppressAutoMessageDuringRefresh = false;
    chromSelect.innerHTML = '<option value="" selected hidden>Error loading chromosomes</option>';
        });
    }

    // When the user opens the chromosome dropdown, show the full species list
    // (replace any earlier "detected-only" list produced after an upload). Use
    // mousedown so we populate options before the native select opens.
    chromSelect.addEventListener('mousedown', function (e) {
      try {
        const species = speciesSelect.value;
        if (!species) return;
        // fetch full list but DO NOT auto-select the server-suggested chromosome
        suppressAutoMessageDuringRefresh = true;
        fetch(`/chromosomes?species=${encodeURIComponent(species)}`)
          .then(r => r.json())
          .then(data => {
            suppressAutoMessageDuringRefresh = false;
            if (data.chromosomes && Array.isArray(data.chromosomes)) {
              populateChromOptions(data.chromosomes, false);
            }
          })
          .catch(() => { suppressAutoMessageDuringRefresh = false; });
      } catch (e) { /* ignore */ }
    });
    speciesSelect.addEventListener('change', updateChromosomes);
  // if user manually chooses a chromosome, forget the server suggestion
  chromSelect.addEventListener('change', () => { try { lastSuggestedChromosome = null; } catch(_) {} });
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

    // If a staged upload exists, submit the analyze form then go straight to processing page
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
  // Explicitly request analysis start for staged uploads (server will not auto-run otherwise)
  fd.append('start_now', '1');
      // POST to /analyze which starts the background job and returns JSON
      fetch(form.action, { method: 'POST', body: fd })
        .then(r => {
          if (!r.ok) throw new Error('Server returned ' + r.status);
          return r.json();
        })
        .then(data => {
          if (!data || !data.run_id) throw new Error('Invalid response from server');
          const next = data.wait_url || (`/runs/${encodeURIComponent(data.run_id)}/final`);
          window.location.href = next;
        })
        .catch(err => {
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