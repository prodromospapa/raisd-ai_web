const dz = document.getElementById('dropzone');
const form = document.getElementById('job-form');
const submitBtn = document.getElementById('submit-btn');

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
    // reevaluate readiness when file selection changes
    try { updateAnalyzeEnabled(); } catch(_) {}
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

// Determine whether the Analyze button should be enabled
function isAnalyzeReady() {
  if (!form) return false;
  const speciesEl = document.getElementById('species-select') || form.querySelector('select[name=species]');
  const chromEl = document.getElementById('chromosome-select') || form.querySelector('select[name=chromosome]');
  const gridEl = form.querySelector('input[name=grid]');
  const species = speciesEl && speciesEl.value;
  const chromosome = chromEl && chromEl.value;
  const gridVal = gridEl && gridEl.value;
  const gridOk = gridVal && Number(gridVal) > 0 && Number.isFinite(Number(gridVal));
  const hasStaged = Boolean(stagedRunId && stagedFilename);
  const fileInput = document.getElementById('file-input');
  const hasLocalFile = !!(fileInput && fileInput.files && fileInput.files.length > 0);
  if (isUploading) return false;
  // If the user selected Combined as the metric, ensure they picked at least two
  try {
    const metricEl = form.querySelector('select[name=distance_metric]') || document.getElementById('distance-metric-select-form');
    const metricVal = metricEl ? String(metricEl.value || '').toLowerCase() : '';
    if (metricVal === 'combined') {
      try {
        // Require the user to have explicitly selected at least two combined metrics.
        // Use the helper which now returns only actually-checked metrics (no fallback).
        const picked = getCheckedCombinedMetricsFromForm();
        try { if (window && window.console && window.console.debug) { console.debug('isAnalyzeReady picked combined metrics:', picked); } }catch(_){ }
        if (!picked || picked.length < 2) return false;
      } catch (e) { /* ignore and proceed */ }
    }
  } catch (e) { /* ignore and proceed */ }
  try { if (window && window.console && window.console.debug) { console.debug('isAnalyzeReady final:', { species: species, chromosome: chromosome, gridOk: gridOk, hasStaged: hasStaged, hasLocalFile: hasLocalFile }); } } catch(_){}
  return Boolean(species && chromosome && gridOk && (hasStaged || hasLocalFile));
}

function updateAnalyzeEnabled() {
  if (!submitBtn) return;
  const ready = isAnalyzeReady();
  submitBtn.disabled = !ready;
  submitBtn.setAttribute('aria-disabled', String(!ready));
}

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
  if (submitBtn) {
    // Disable during upload or when inputs are incomplete
    submitBtn.disabled = isUploading || !isAnalyzeReady();
    submitBtn.setAttribute('aria-disabled', String(submitBtn.disabled));
  }
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
  try { updateAnalyzeEnabled(); } catch(_) {}
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
    try { updateAnalyzeEnabled(); } catch(_) {}
  };
  xhr.onload = function () {
  isUploading = false;
    lastUploadXhr = null;
    try {
      let data = xhr.response;
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
        try { updateAnalyzeEnabled(); } catch(_) {}
      } else if (data && data.error) {
  alert('Upload failed: ' + data.error);
        resetDropzone();
        try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
        try { updateAnalyzeEnabled(); } catch(_) {}
      } else if (xhr.status < 200 || xhr.status >= 300) {
        alert('Upload failed: server returned ' + xhr.status);
        resetDropzone();
        try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
        try { updateAnalyzeEnabled(); } catch(_) {}
      }

      // ------------------ All Distances / Metric selector UI ------------------
      // Moved from inline template for better separation and maintainability.
      window.initAllDistancesUI = function initAllDistancesUI(){
        try{
          const md = JSON.parse(document.getElementById('match-data')?.textContent || '{}') || {};
          const allDistances = (md.all_distances && Array.isArray(md.all_distances)) ? md.all_distances.map(x=>({demographic_model:x.demographic_model,population:x.population,distances:x.distances})) : null;
          let metricsAvailable = md.metrics_available || (allDistances && allDistances.length ? Object.keys(allDistances[0].distances || {}) : ['JSD']);
          // sort metricsAvailable alphanumerically case-insensitive
          try{ metricsAvailable = metricsAvailable.slice().sort((a,b)=> String(a).toLowerCase().localeCompare(String(b).toLowerCase())); }catch(_){ }
          // Ensure the server-selected metric is present in the options so the UI
          // can pick it even when metrics_available was not present or incomplete.
          try {
            const sel = (md && md.selected_metric) ? md.selected_metric : (md && md.match && md.match.selected_metric ? md.match.selected_metric : null);
            if (sel && metricsAvailable && Array.isArray(metricsAvailable) && !metricsAvailable.includes(sel)) {
              metricsAvailable = [sel].concat(metricsAvailable);
            }
          } catch (e) { /* ignore */ }
          const metricSelect = document.getElementById('distance-metric-select');
          const combinedSelect = document.getElementById('combined-metric-select');
          const combinedCheckboxes = document.getElementById('combined-metric-checkboxes');
          const tableBody = document.querySelector('#allMatchesTable tbody');
          if(!metricSelect || !tableBody) return;
          // populate options - avoid clobbering server-rendered options when possible.
          try {
            // Always ensure Combined is first and the remaining metrics
            // are sorted case-insensitively. Preserve current selection if present.
            try{
              const toShow = (metricsAvailable || ['JSD']).slice().filter(v=>v);
              // dedupe and sort
              const uniq = Array.from(new Set(toShow)).sort((a,b)=> String(a).toLowerCase().localeCompare(String(b).toLowerCase()));
              const finalList = ['Combined'].concat(uniq.filter(x=>String(x) !== 'Combined'));
              const curVal = metricSelect.value;
              metricSelect.innerHTML = '';
              finalList.forEach(m=>{ const o=document.createElement('option'); o.value=m; o.text=m; metricSelect.appendChild(o); });
              if(curVal && Array.from(metricSelect.options).some(o=>o.value===curVal)) metricSelect.value = curVal;
            }catch(_){ }
          } catch (e) { /* ignore */ }
          // Default selection: prefer the metric used for the best match if provided by server
          try {
            // Re-read merged match-data to prefer any values provided by match_details
            const md_local = JSON.parse(document.getElementById('match-data')?.textContent || '{}') || {};
            const selMetric = (md_local && md_local.selected_metric) ? md_local.selected_metric : (md_local && md_local.match && md_local.match.selected_metric ? md_local.match.selected_metric : null);
            if (selMetric) {
              // If the metric exists in the options, pick it; else default to 'JSD'
              const found = Array.from(metricSelect.options).some(o => o.value === selMetric);
              if (found) metricSelect.value = selMetric;
              else metricSelect.value = 'JSD';
            }
          } catch (e) { /* ignore */ }
          if(combinedSelect){ combinedSelect.innerHTML=''; // keep Combined out of combinedSelect options; it represents the multi-metric picker
            (metricsAvailable||[]).forEach(m=>{ const o=document.createElement('option'); o.value=m; o.text=m; combinedSelect.appendChild(o); }); }
          if (combinedCheckboxes) {
            combinedCheckboxes.innerHTML = '';
            // Controls row: Select all / Clear and count
            const controls = document.createElement('div'); controls.className = 'combined-controls';
            const selectAll = document.createElement('button'); selectAll.type='button'; selectAll.className='combined-control-btn'; selectAll.textContent = 'Select all';
            const clearBtn = document.createElement('button'); clearBtn.type='button'; clearBtn.className='combined-control-btn'; clearBtn.textContent = 'Clear';
            const countSpan = document.createElement('div'); countSpan.className = 'combined-count'; countSpan.textContent = (metricsAvailable||[]).length + ' selected';
            controls.appendChild(selectAll); controls.appendChild(clearBtn); controls.appendChild(countSpan);
            combinedCheckboxes.appendChild(controls);

            const pills = document.createElement('div'); pills.className = 'combined-pills';
            const state = { selected: new Set(metricsAvailable || []) };
            (metricsAvailable||[]).forEach(m => {
              const id = 'combined_pill_' + m.replace(/[^a-z0-9]/gi,'_');
              const btn = document.createElement('button'); btn.type='button'; btn.className='combined-pill active'; btn.setAttribute('data-metric', m); btn.id = id; btn.title = m;
              const chk = document.createElement('input'); chk.type='checkbox'; chk.value = m; chk.checked = true;
              const span = document.createElement('span'); span.textContent = m;
              btn.appendChild(chk); btn.appendChild(span);
              btn.addEventListener('click', ()=>{
                const val = btn.getAttribute('data-metric');
                if(state.selected.has(val)){ state.selected.delete(val); btn.classList.remove('active'); chk.checked = false; }
                else { state.selected.add(val); btn.classList.add('active'); chk.checked = true; }
                countSpan.textContent = state.selected.size + ' selected';
                // fire a custom change event on the container so other listeners update
                try{ combinedCheckboxes.dispatchEvent(new Event('change')); }catch(_){ }
              });
              pills.appendChild(btn);
            });
            combinedCheckboxes.appendChild(pills);

            selectAll.addEventListener('click', ()=>{
              (metricsAvailable||[]).forEach(m=> state.selected.add(m));
              Array.from(pills.children).forEach(b=>{ b.classList.add('active'); b.querySelector('input').checked = true; });
              countSpan.textContent = state.selected.size + ' selected';
              try{ combinedCheckboxes.dispatchEvent(new Event('change')); }catch(_){ }
            });
            clearBtn.addEventListener('click', ()=>{
              state.selected.clear();
              Array.from(pills.children).forEach(b=>{ b.classList.remove('active'); b.querySelector('input').checked = false; });
              countSpan.textContent = '0 selected';
              try{ combinedCheckboxes.dispatchEvent(new Event('change')); }catch(_){ }
            });
            // ensure there's a toggle button in the DOM (for template this may be pre-rendered)
            try {
              const runMeta = document.getElementById('run-meta');
              const runId = runMeta ? runMeta.getAttribute('data-run-id') : null;
              const toggleBtn = document.getElementById('combined-toggle-btn');
              const storageKey = runId ? ('combined_visible_' + runId) : 'combined_visible';
              // helper to set visibility and persist. Use the outer container as authoritative
              function setCombinedVisible(vis, persist=true){
                try{ const outer = document.getElementById('combined-metric-container'); if(outer) outer.style.display = vis ? '' : 'none'; }catch(_){ }
                try{ if(combinedCheckboxes) combinedCheckboxes.style.display = vis ? 'block' : 'none'; }catch(_){}
                if(toggleBtn){ toggleBtn.textContent = vis ? 'Hide' : 'Show'; toggleBtn.setAttribute('aria-expanded', String(vis)); }
                try{ if(persist){ localStorage.setItem(storageKey, vis ? '1' : '0'); } }catch(_){ }
                try{ if(typeof adjustAllMatchesScroll === 'function') adjustAllMatchesScroll(); }catch(_){ }
              }
              // Expose globally so other inline handlers can call it
              try{ window.setCombinedVisible = setCombinedVisible; }catch(_){ }
              // initialize from storage (default: visible when Combined selected)
              try{
                const metricSel = document.getElementById('distance-metric-select');
                const isCombined = metricSel && String(metricSel.value).toLowerCase() === 'combined';
                const stored = (function(){ try{ return localStorage.getItem(storageKey); }catch(_){ return null; }})();
                if(isCombined && stored !== null){ setCombinedVisible(stored === '1', false); }
                else { // no stored pref or not Combined -> reflect current container/checkbox state
                  const outer = document.getElementById('combined-metric-container');
                  const cur = outer ? (outer.style.display !== 'none') : true;
                  setCombinedVisible(cur, false);
                }
              }catch(_){ setCombinedVisible(true, false); }
              if(toggleBtn){ toggleBtn.addEventListener('click', ()=>{ try{ const outer = document.getElementById('combined-metric-container'); const cur = outer ? (outer.style.display !== 'none') : false; setCombinedVisible(!cur, true); }catch(_){ } }); }
            } catch (e) { /* ignore toggle wiring errors */ }
          }
          // default: select all metrics for Combined
          if (combinedSelect) {
            try {
              Array.from(combinedSelect.options).forEach(o => { o.selected = true; });
            } catch (e) { /* ignore */ }
          }

          function computeCombined(list, metrics){
            const N = list.length; if(!N) return [];
            const ranks = {};
            metrics.forEach(m=>{ ranks[m]=new Array(N); const arr = list.map((it,idx)=>({idx:idx,val:(it.distances&&it.distances[m]!=null)?it.distances[m]:Infinity})); arr.sort((a,b)=> (a.val===b.val)? a.idx-b.idx : a.val - b.val); for(let i=0;i<arr.length;i++){ ranks[m][arr[i].idx]=i+1; } });
            const avg = new Array(N); for(let i=0;i<N;i++){ let s=0; metrics.forEach(m=>{ s += (ranks[m][i]||N); }); avg[i]=s/metrics.length; }
            return avg.map(v=>v/N);
          }

          function renderByMetric(metric, combinedMetrics){
            // Re-read match-data each time so the table always reflects the
            // most up-to-date merged payload (including arrays from match_details).
            const md_local = JSON.parse(document.getElementById('match-data')?.textContent || '{}') || {};
            // Prefer all_distances (per-metric distances) when available; fall back to all_jsd.
            const source = (md_local.all_distances && Array.isArray(md_local.all_distances)) ? md_local.all_distances : (md_local.all_jsd && Array.isArray(md_local.all_jsd) ? md_local.all_jsd : []);
            const list = source.map(x => (x.distances ? x : {demographic_model: x.demographic_model, population: x.population, distances: {JSD: x.jsd}}));
            // Ensure metricsAvailable is derived from current data when possible
            const localMetricsAvailable = md_local.metrics_available || (list.length ? Object.keys(list[0].distances || {}) : ['JSD']);
            let order = list.map((_,i)=>i);
            if(metric === 'Combined'){
              const metrics = (combinedMetrics && combinedMetrics.length) ? combinedMetrics : (localMetricsAvailable.slice ? localMetricsAvailable.slice() : []);
              // If no metrics selected for Combined, clear table and return
              if(!metrics || metrics.length === 0){
                try { const hdr = document.getElementById('allMatchesMetricHeader'); if(hdr) hdr.textContent = 'Combined'; } catch(_) {}
                tableBody.innerHTML = '';
                return;
              }
              const norm = computeCombined(list, metrics);
              order.sort((a,b)=> norm[a] - norm[b]);
            } else {
              order.sort((a,b)=> { const va=(list[a].distances&&list[a].distances[metric]!=null)?list[a].distances[metric]:Infinity; const vb=(list[b].distances&&list[b].distances[metric]!=null)?list[b].distances[metric]:Infinity; if(va===vb) return a-b; return va-vb; });
            }
            // debug: report which metric we are rendering and whether distances exist
            try {
              if (window && window.console && window.console.log) {
                const sample = list && list.length ? list[0] : null;
                const hasMetric = sample && sample.distances && (sample.distances[metric] != null);
                console.log('renderByMetric:', { metric: metric, hasMetric: !!hasMetric, combinedMetrics: combinedMetrics, sampleMetricKeys: sample?Object.keys(sample.distances||{}):[] });
              }
            } catch (e) { /* noop */ }
            // update header label for the metric column
            try {
              const hdr = document.getElementById('allMatchesMetricHeader');
              if (hdr) hdr.textContent = metric === 'Combined' ? 'Combined' : metric;
            } catch (e) { /* ignore */ }
            tableBody.innerHTML='';
            order.forEach((idx,i)=>{
              const m = list[idx];
              const tr = document.createElement('tr'); tr.setAttribute('data-dm', m.demographic_model||''); tr.setAttribute('data-pop', m.population||''); tr.style.cursor='pointer';
              const idxTd = document.createElement('td'); idxTd.style.padding='6px'; idxTd.textContent = String(i+1);
              const mod = document.createElement('td'); mod.style.padding='6px'; mod.textContent = m.demographic_model||'';
              const pop = document.createElement('td'); pop.style.padding='6px'; pop.textContent = m.population||'';
              const valTd = document.createElement('td'); valTd.style.padding='6px'; valTd.style.textAlign='right';
              let display = '';
              if(metric === 'Combined'){
                const metrics = (combinedMetrics && combinedMetrics.length)? combinedMetrics : (localMetricsAvailable.slice ? localMetricsAvailable.slice() : []);
                const norm = computeCombined(list, metrics);
                display = isFinite(norm[idx])? norm[idx].toFixed(6) : '';
              } else {
                const v = (m.distances && m.distances[metric] != null)? m.distances[metric] : null;
                display = (v !== null && v !== undefined && isFinite(v)) ? Number(v).toFixed(6) : '';
              }
              valTd.textContent = display;
              tr.appendChild(idxTd); tr.appendChild(mod); tr.appendChild(pop); tr.appendChild(valTd);
              tr.addEventListener('click', ()=> startRescan(tr, m.demographic_model, m.population));
              tr.addEventListener('keypress', (e)=>{ if(e.key==='Enter' || e.key===' '){ e.preventDefault(); startRescan(tr, m.demographic_model, m.population); }});
              tr.setAttribute('tabindex','0');
              tableBody.appendChild(tr);
            });
          }

          metricSelect.addEventListener('change', ()=>{
            const v = metricSelect.value;
            // update header immediately so UI reflects the selected metric without waiting
            try {
              const hdr = document.getElementById('allMatchesMetricHeader');
              if (hdr) hdr.textContent = v === 'Combined' ? 'Combined' : v;
            } catch (e) { /* ignore */ }
            if(v === 'Combined'){
              // When switching to Combined, always expose both the toggle button and
              // the combined metrics list so the user can immediately see and edit
              // which metrics are included. Use the helper setCombinedVisible if
              // available (it will avoid persisting this forced visibility).
              try{
                if(typeof setCombinedVisible === 'function'){
                  try{ setCombinedVisible(true, false); }catch(_){ /* ignore */ }
                } else {
                  if(combinedCheckboxes) combinedCheckboxes.style.display = 'block';
                  if(toggleBtn){ toggleBtn.textContent = 'Hide'; toggleBtn.setAttribute('aria-expanded', 'true'); }
                }
              }catch(_){ if(combinedCheckboxes) combinedCheckboxes.style.display = 'block'; }
              combinedSelect.style.display = 'inline-block';
              // ensure all options are selected by default when Combined is chosen
              try { Array.from(combinedSelect.options).forEach(o => { o.selected = true; }); } catch(e){}
            } else if(combinedSelect){ combinedSelect.style.display = 'none'; }
            const sel = (combinedSelect? Array.from(combinedSelect.options).filter(o=>o.selected).map(o=>o.value) : []);
            // call render to compute values for the chosen metric
            renderByMetric(v, sel);
          });
          if (combinedCheckboxes) {
            // Listen for delegated change events (pills dispatch a change on toggle)
            combinedCheckboxes.addEventListener('change', ()=>{
              if(metricSelect.value === 'Combined'){
                const selElems = Array.from(combinedCheckboxes.querySelectorAll('input[type=checkbox]'));
                const sel = selElems.filter(i=>i.checked).map(i=>i.value);
                // Persistent feedback: show message until user picks >=2 metrics
                try{
                  const msg = combinedCheckboxes.querySelector('.combined-feedback');
                  if(msg){ msg.textContent = (sel.length < 2) ? 'Pick at least two metrics to enable Combined' : ''; }
                }catch(_){ }
                renderByMetric('Combined', sel);
              }
            });
            // Also handle keyboard activation on pill buttons (space/enter)
            combinedCheckboxes.addEventListener('keydown', (e)=>{
              if(e.key === 'Enter' || e.key === ' '){ const btn = document.activeElement; if(btn && btn.classList && btn.classList.contains('combined-pill')){ btn.click(); e.preventDefault(); } }
            });
          }
          // initial render: default to first metric or JSD
          // initial: choose metric and default Combined selection
          const initMetric = metricSelect.value || 'JSD';
          if (initMetric === 'Combined') {
            // ensure all checkboxes are checked
            if (combinedCheckboxes) {
              Array.from(combinedCheckboxes.querySelectorAll('input[type=checkbox]')).forEach(i=>{ i.checked = true; });
              const sel = Array.from(combinedCheckboxes.querySelectorAll('input[type=checkbox]')).filter(i=>i.checked).map(i=>i.value);
              // initial feedback: if less than two selected, show persistent message
              try{ const msg = combinedCheckboxes.querySelector('.combined-feedback'); if(msg){ msg.textContent = (sel.length < 2) ? 'Pick at least two metrics to enable Combined' : ''; } }catch(_){ }
              renderByMetric('Combined', sel);
            } else if (combinedSelect) {
              try { Array.from(combinedSelect.options).forEach(o => { o.selected = true; }); } catch(e){}
              const sel = Array.from(combinedSelect.options).filter(o=>o.selected).map(o=>o.value);
              renderByMetric('Combined', sel);
            } else {
              renderByMetric('Combined', []);
            }
          } else {
            renderByMetric(initMetric, []);
          }
        }catch(e){ console.warn('metric selector init failed', e); }
      };

          // Auto-init on pages that include match-data
          try{
            if(document.getElementById && document.getElementById('match-data')){
              window.addEventListener('DOMContentLoaded', async ()=>{
                try{
                  // If match-data is missing heavy arrays, try fetching them lazily
                  const mdEl = document.getElementById('match-data');
                  let md = {};
                  try { md = JSON.parse(mdEl.textContent || '{}') || {}; } catch(_) { md = {}; }
                  const runMeta = document.getElementById('run-meta');
                  const runId = runMeta ? runMeta.getAttribute('data-run-id') : null;
                  let needHeavy = false;
                  if(!md || (typeof md !== 'object')) md = {};
                  if((!md.all_distances || !Array.isArray(md.all_distances)) || (!md.input_sfs && !md.observed_sfs)){
                    needHeavy = Boolean(runId);
                  }
                  if(needHeavy && runId){
                    try{
                      const resp = await fetch(`/runs/${encodeURIComponent(runId)}/match_heavy`);
                      if(resp && resp.ok){
                        // parse JSON once (browser handles gzip/encoding)
                        let payload = null;
                        try { payload = await resp.json(); } catch(err){ payload = null; }
                        if(payload){
                          // Merge payload.match into md first
                          if(payload.match && typeof payload.match === 'object'){
                              md = Object.assign({}, md, payload.match);
                          }
                          // Merge match_details (heavy arrays) into md when present
                          if(payload.match_details && typeof payload.match_details === 'object'){
                            const d = payload.match_details;
                              // Prefer arrays from match_details (they are the full, un-pruned versions).
                              if(Array.isArray(d.observed_sfs)) md.input_sfs = d.observed_sfs;
                              if(Array.isArray(d.best_expected_sfs)) md.best_expected_sfs = d.best_expected_sfs;
                              if(Array.isArray(d.top_matches)) md.top_matches = d.top_matches;
                              if(Array.isArray(d.all_jsd)) md.all_jsd = d.all_jsd;
                              if(Array.isArray(d.all_distances)) md.all_distances = d.all_distances;
                              if(d.metrics_available && Array.isArray(d.metrics_available)) md.metrics_available = d.metrics_available;
                          }
                          // Persist merged object back into the DOM so initializers read updated data
                          try { mdEl.textContent = JSON.stringify(md); } catch(e){ /* ignore */ }
                          // Debug: expose merged metrics info to the console for troubleshooting
                          try {
                            if (window && window.console && window.console.log) {
                              const info = { metrics_available: md.metrics_available, selected_metric: md.selected_metric, has_all_distances: Array.isArray(md.all_distances) };
                              console.log('match_heavy merged:', info);
                            }
                          } catch (e) { /* noop */ }
                        }
                      }
                    }catch(e){ /* best-effort; continue */ }
                  }
                  // Now initialize AllDistances UI and hover tooltips (they will re-read match-data)
                  try{ if(window.initAllDistancesUI) window.initAllDistancesUI(); }catch(e){}
                  try{ if(window.initHoverTooltips) window.initHoverTooltips(); }catch(e){}
                }catch(_){ }
              });
            }
          }catch(_){ }
    } catch (e) {
      alert('Upload failed: ' + e);
      resetDropzone();
      try { const bar = document.getElementById('upload-bar'); if (bar) bar.classList.remove('animate'); } catch(_) {}
      try { updateAnalyzeEnabled(); } catch(_) {}
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
  // include selected metric preferences
  try {
    const metric = form.querySelector('select[name=distance_metric]');
    if (metric && metric.value) data.append('distance_metric', metric.value);
    const combinedValues = getCheckedCombinedMetricsFromForm();
    combinedValues.forEach(v => data.append('combined_metrics', v));
  } catch (e) { /* ignore */ }

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
    try { updateAnalyzeEnabled(); } catch(_) {}
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
  try { updateAnalyzeEnabled(); } catch(_) {}

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
        try { updateAnalyzeEnabled(); } catch(_) {}
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
          try { updateAnalyzeEnabled(); } catch(_) {}
        })
        .catch(() => {
          suppressAutoMessageDuringRefresh = false;
    chromSelect.innerHTML = '<option value="" selected hidden>Error loading chromosomes</option>';
          try { updateAnalyzeEnabled(); } catch(_) {}
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
            try { updateAnalyzeEnabled(); } catch(_) {}
          })
          .catch(() => { suppressAutoMessageDuringRefresh = false; try { updateAnalyzeEnabled(); } catch(_) {} });
      } catch (e) { /* ignore */ }
    });
    speciesSelect.addEventListener('change', () => { updateChromosomes(); try { updateAnalyzeEnabled(); } catch(_) {} });
  // if user manually chooses a chromosome, forget the server suggestion
  chromSelect.addEventListener('change', () => { try { lastSuggestedChromosome = null; } catch(_) {}; try { updateAnalyzeEnabled(); } catch(_) {} });
    // Initial population on page load
    updateChromosomes();
  }

  // Show/hide combined-metrics control on the index page based on selected distance metric
  try {
    const metricFormSelect = document.getElementById('distance-metric-select-form');
    const combinedWrapper = document.getElementById('combined-metrics-wrapper');
    if (metricFormSelect && combinedWrapper) {
      function populateCombinedCheckboxes(containerEl) {
        try {
          // Prefer match-data.metrics_available (keeps server-provided order) like the results page.
          let metrics = [];
          let metricsFromMatchData = false;
          try {
            const md = JSON.parse(document.getElementById('match-data')?.textContent || '{}') || {};
            if (md.metrics_available && Array.isArray(md.metrics_available) && md.metrics_available.length) {
              metrics = md.metrics_available.slice();
              metricsFromMatchData = true;
            }
          } catch(_) { /* ignore parse errors */ }
          if (!metricsFromMatchData) {
            metrics = Array.from(document.querySelectorAll('select[name=distance_metric] option')).map(o=>o.value).filter(v=>v && v !== 'Combined');
            // sort available metrics for display when derived locally
            try{ metrics = metrics.slice().sort((a,b)=> String(a).toLowerCase().localeCompare(String(b).toLowerCase())); }catch(_){ }
          }
          if(!containerEl) return metrics;
          containerEl.innerHTML = '';
          // Top controls: Select All / Clear and count (mirrors result.html)
          const controls = document.createElement('div'); controls.className = 'combined-controls';
          const selectAll = document.createElement('button'); selectAll.type='button'; selectAll.className='combined-control-btn'; selectAll.textContent='Select all';
          const clearBtn = document.createElement('button'); clearBtn.type='button'; clearBtn.className='combined-control-btn'; clearBtn.textContent='Clear';
          const countSpan = document.createElement('div'); countSpan.className = 'combined-count'; countSpan.textContent = metrics.length + ' selected';
          controls.appendChild(selectAll); controls.appendChild(clearBtn); controls.appendChild(countSpan);
          containerEl.appendChild(controls);

          const pills = document.createElement('div'); pills.className = 'combined-pills';
          const state = { selected: new Set(metrics) };
          // feedback message
          const feedback = document.createElement('div'); feedback.className = 'combined-feedback'; feedback.style.fontSize='12px'; feedback.style.color='var(--muted)'; feedback.style.padding='6px 0 0 0';
          try{ if(state.selected.size < 2){ feedback.textContent = 'Pick at least two metrics to enable Combined'; } }catch(_){ }

          // create pills (all active by default)
          metrics.forEach(m => {
            const id = 'combined_pill_' + m.replace(/[^a-z0-9]/gi,'_');
            const btn = document.createElement('button'); btn.type='button'; btn.className='combined-pill active'; btn.setAttribute('data-metric', m); btn.id = id; btn.title = m;
            const chk = document.createElement('input'); chk.type='checkbox'; chk.value = m; chk.name = 'combined_metrics'; chk.checked = true; chk.style.display = 'none';
            const span = document.createElement('span'); span.textContent = m;
            btn.appendChild(chk); btn.appendChild(span);
            btn.addEventListener('click', ()=>{
              const val = btn.getAttribute('data-metric');
              if(state.selected.has(val)){ state.selected.delete(val); btn.classList.remove('active'); chk.checked = false; }
              else { state.selected.add(val); btn.classList.add('active'); chk.checked = true; }
              countSpan.textContent = state.selected.size + ' selected';
              try{ containerEl.dispatchEvent(new Event('change')); }catch(_){ }
              try{ updateAnalyzeEnabled(); }catch(_){}
              try{ feedback.textContent = (state.selected.size < 2) ? 'Pick at least two metrics to enable Combined' : ''; }catch(_){ }
            });
            pills.appendChild(btn);
          });
          containerEl.appendChild(pills);
          containerEl.appendChild(feedback);
          try{ updateAnalyzeEnabled(); }catch(_){}

          // wire control buttons
          selectAll.addEventListener('click', (e)=>{ e.preventDefault(); metrics.forEach(m=> state.selected.add(m)); Array.from(pills.children).forEach(b=>{ b.classList.add('active'); b.querySelector('input').checked = true; }); countSpan.textContent = state.selected.size + ' selected'; try{ containerEl.dispatchEvent(new Event('change')); }catch(_){ } try{ feedback.textContent = ''; }catch(_){ } try{ updateAnalyzeEnabled(); }catch(_){} });
          clearBtn.addEventListener('click', (e)=>{ e.preventDefault(); state.selected.clear(); Array.from(pills.children).forEach(b=>{ b.classList.remove('active'); b.querySelector('input').checked = false; }); countSpan.textContent = '0 selected'; try{ containerEl.dispatchEvent(new Event('change')); }catch(_){ } try{ feedback.textContent = 'Pick at least two metrics to enable Combined'; }catch(_){ } try{ updateAnalyzeEnabled(); }catch(_){} });

          return metrics;
        } catch (e) { /* ignore */ return []; }
      }

      function updateCombinedVisibility() {
        try {
          const isCombined = metricFormSelect.value === 'Combined';
          if (isCombined) {
            combinedWrapper.style.display = 'block';
            // shrink the metric select so the toggle fits to its right
            try { metricFormSelect.classList.add('compact-metric-select'); } catch(_){}
            // show the toggle button
            try { const toggleBtn = document.getElementById('combined-toggle-btn'); if(toggleBtn) toggleBtn.style.display = ''; } catch(_){}
            // Prefer updating the dedicated combined checkbox container when present
            const cbBox = document.getElementById('combined-metric-checkboxes');
            const legacy = document.getElementById('combined-metrics-checkboxes-form');
            if (cbBox) {
              cbBox.style.display = 'block';
              populateCombinedCheckboxes(cbBox);
            } else if (legacy) {
              legacy.style.display = 'flex';
              populateCombinedCheckboxes(legacy);
            }
            try{ updateAnalyzeEnabled(); }catch(_){}
          } else {
            combinedWrapper.style.display = 'none';
            // restore select width and hide toggle
            try { metricFormSelect.classList.remove('compact-metric-select'); } catch(_){}
            try { const toggleBtn = document.getElementById('combined-toggle-btn'); if(toggleBtn) toggleBtn.style.display = 'none'; } catch(_){}
            const cbBox = document.getElementById('combined-metric-checkboxes'); if (cbBox) cbBox.style.display = 'none';
            const legacy = document.getElementById('combined-metrics-checkboxes-form'); if (legacy) legacy.style.display = 'none';
          }
        } catch (e) { /* ignore */ }
      }
      metricFormSelect.addEventListener('change', updateCombinedVisibility);
      // initial
      updateCombinedVisibility();
    }
  } catch (e) { /* ignore */ }

  // Helper: read checked combined metrics from either checkbox container or legacy multi-select
  function getCheckedCombinedMetricsFromForm() {
    // Prefer the canonical combined checkbox container (used on results/index)
    let cbContainer = document.getElementById('combined-metric-checkboxes');
    if (!cbContainer) cbContainer = document.getElementById('combined-metrics-checkboxes-form');
    if (cbContainer && cbContainer.querySelectorAll) {
      const checks = Array.from(cbContainer.querySelectorAll('input[type=checkbox][name="combined_metrics"]'));
      if (checks && checks.length){
        const picked = checks.filter(c=>c.checked).map(c=>c.value);
        try{ if(window && window.console && window.console.debug) console.debug('getCheckedCombinedMetricsFromForm: found checkboxes, picked=', picked); }catch(_){ }
        return picked;
      }
      // If the container exists but no checkboxes rendered yet, treat as no selection (require explicit picks)
      try{ if(window && window.console && window.console.debug) console.debug('getCheckedCombinedMetricsFromForm: container exists but no checkboxes rendered yet'); }catch(_){ }
      return [];
    }
    // fallback: legacy multi-select
    const combinedSel = document.querySelector('select[name=combined_metrics]');
    if (combinedSel) return Array.from(combinedSel.options).filter(o=>o.selected).map(o=>o.value);
    return [];
  }
}

// Keep submit behavior: if user presses Analyze and a file is present, start upload (or ignore if already uploading)
if (form) {
  // watch grid validity to toggle readiness
  try {
    const gridEl = form.querySelector('input[name=grid]');
    if (gridEl) ['input','change','blur'].forEach(ev => gridEl.addEventListener(ev, () => { try { updateAnalyzeEnabled(); } catch(_) {} }));
  } catch(_) {}
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
      // include user's selected metric preferences so the background worker uses them
      try {
        const metric = form.querySelector('select[name=distance_metric]');
        if (metric && metric.value) fd.append('distance_metric', metric.value);
        const combinedValues = getCheckedCombinedMetricsFromForm();
        if(metric && metric.value === 'Combined' && combinedValues.length < 2){ alert('Please select at least two metrics for Combined.'); return; }
        combinedValues.forEach(v => fd.append('combined_metrics', v));
      } catch (e) { /* ignore */ }
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
  // Validate Combined selection on client-side before upload
  try{
    const metric = form.querySelector('select[name=distance_metric]');
    if(metric && metric.value === 'Combined'){
      const combinedValues = getCheckedCombinedMetricsFromForm();
      if(!combinedValues || combinedValues.length < 2){
        alert('Please select at least two metrics for Combined before analyzing.');
        return;
      }
    }
  }catch(_){ }
  if (!isUploading) uploadFile(file);
  });
}

// Initial state
try { updateAnalyzeEnabled(); } catch(_) {}

// Ensure the index page's combined toggle button works: show/hide the combined metrics
window.addEventListener('DOMContentLoaded', () => {
  try {
    const toggleBtn = document.getElementById('combined-toggle-btn');
    if (!toggleBtn) return;
    const storageKey = 'combined_visible';
    function setCombinedVisibleLocal(vis, persist = true) {
      try { const outer = document.getElementById('combined-metric-container'); if (outer) outer.style.display = vis ? '' : 'none'; } catch(_){}
      try { const cb = document.getElementById('combined-metric-checkboxes'); if (cb) cb.style.display = vis ? 'block' : 'none'; } catch(_){}
      try { toggleBtn.textContent = vis ? 'Hide' : 'Show'; toggleBtn.setAttribute('aria-expanded', String(vis)); } catch(_){}
      try { if (persist) localStorage.setItem(storageKey, vis ? '1' : '0'); } catch(_){}
      try { updateAnalyzeEnabled(); } catch(_){}
    }
    toggleBtn.addEventListener('click', () => {
      try {
        const outer = document.getElementById('combined-metric-container');
        const cur = outer ? (outer.style.display !== 'none') : false;
        setCombinedVisibleLocal(!cur, true);
      } catch(_) {}
    });
    try { const stored = localStorage.getItem(storageKey); if (stored !== null) setCombinedVisibleLocal(stored === '1', false); } catch(_){}
  } catch(_) {}
});