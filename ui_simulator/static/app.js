
// helpers
const $ = sel => document.querySelector(sel);
const $$ = sel => Array.from(document.querySelectorAll(sel));

let ST = window.__INITIAL_STATE__ || {};
let DV = window.__DERIVED__ || {};
let ARTS = {}; // live artifact map for mid-run discovery

function buildCatalogURL(species_id, model_id){
  if(!species_id || !model_id) return null;
  const s = species_id.toLowerCase().replace(/\//g, "_");
  const m = model_id.toLowerCase().replace(/[^a-z0-9_]+/g, "_");
  return `https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_${s}_models_${m}`;
}

async function apiSet(updates){
  const res = await fetch("/api/set",{method:"POST", headers:{"Content-Type":"application/json"}, body: JSON.stringify(updates)});
  const js = await res.json();
  if(!js.ok) throw new Error(js.error || "Update failed");
  ST = js.stt; DV = js.derived;
  renderAll(); // native refresh without page reload
}

function setEnable(el, on){
  if(!el) return;
  if(on){ el.removeAttribute("disabled"); }
  else{ el.setAttribute("disabled",""); }
}
function show(el, on){ if(!el) return; el.style.display = on ? "" : "none"; }
function toggleHide(el, hide){ if(!el) return; el.classList.toggle("hide", !!hide); }

function renderSpeciesModelsContigs(){
  // species list
  const spSel = $("#species");
  spSel.innerHTML = `<option value="">-- select species --</option>` + DV.species_items.map(s=>`<option value="${s.id}" ${ST.species_id===s.id?'selected':''}>${s.label}</option>`).join("");
  // model list
  const mSel = $("#model");
  mSel.innerHTML = `<option value="">-- select model --</option>` + (DV.model_options||[]).map(m=>`<option ${ST.model_id===m?'selected':''}>${m}</option>`).join("");
  setEnable(mSel, !!ST.species_id);
  // contigs
  const cSel = $("#chromosome");
  cSel.innerHTML = `<option value="">-- select chromosome --</option>` + (DV.contigs||[]).map(c=>`<option ${ST.chromosome===c?'selected':''}>${c}</option>`).join("");
  setEnable(cSel, !!(ST.species_id && ST.model_id));
  // catalog link
  const link = $("#catalogLink");
  const url = buildCatalogURL(ST.species_id, ST.model_id);
  if(url){ link.href = url; show(link, true); } else { show(link, false); }
  // chrom len hint
  const hint = $("#chromLenHint");
  if(ST.chromosome && DV.contig_len){
    hint.textContent = `Chromosome ${ST.chromosome} length: ${Number(DV.contig_len).toLocaleString()} bp`;
  } else { hint.textContent = ""; }
}

function renderPops(){
  const box = $("#popBuilder");
  const pops = ST.populations || [];
  if(!(ST.species_id && ST.model_id && ST.chromosome)){
    box.innerHTML = `<div class="hint">Select previous steps to configure populations.</div>`;
    return;
  }
  box.innerHTML = "";
  pops.forEach((p,i)=>{
    const row = document.createElement("div");
    row.className = "pop-row";
    row.innerHTML = `
      <label class="toggle"><input type="checkbox" ${p.selected?'checked':''} data-key="pop_selected_${i}"> ${p.name}</label>
      <input type="number" min="0" value="${p.n||0}" data-key="pop_n_${i}">
    `;
    box.appendChild(row);
  });
  box.querySelectorAll("input[type=checkbox]").forEach(chk=>{
    chk.addEventListener("change", e=> apiSet({[e.target.dataset.key]: e.target.checked}));
  });
  box.querySelectorAll("input[type=number]").forEach(nm=>{
    nm.addEventListener("change", e=> apiSet({[e.target.dataset.key]: e.target.value}));
  });
}

function renderCore(){
  const ready = DV.ordered_ready;
  // Render engine options dynamically so we can mark the currently selected
  // sweep engine with an asterisk that moves when the user changes engine.
  const engineSel = $("#engine");
  const engines = ['discoal','ms','msms','scrm','msprime'];
  // Render engine options without the asterisk in the core dropdown
  engineSel.innerHTML = engines.map(e => `<option value="${e}" ${ST.engine===e?'selected':''}>${e}</option>`).join('');
  setEnable(engineSel, ready);
  setEnable($("#replicates"), ready);

  const chromNative = $("#chrom_length_mode");
  chromNative.checked = !!ST.chrom_length_mode;
  setEnable(chromNative, ready);

  const seqLen = $("#seq_length");
  const tgtSnps = $("#target_snps");
  const tolAuto = $("#tol_auto");
  const tolVal = $("#target_snps_tol");

  seqLen.value = ST.seq_length||0;
  tgtSnps.value = ST.target_snps||0;
  tolAuto.checked = !!ST.target_snps_tol_auto;
  tolVal.value = ST.target_snps_tol||10;

  setEnable(seqLen, ready && !ST.chrom_length_mode && (ST.target_snps<=0));
  setEnable(tgtSnps, ready && !ST.chrom_length_mode && (ST.seq_length<=0));
  setEnable(tolAuto, ready && !ST.chrom_length_mode && (ST.target_snps>0));
  setEnable(tolVal, ready && !ST.chrom_length_mode && (ST.target_snps>0) && !ST.target_snps_tol_auto);

  // output
  $("#output_format").value = ST.output_format||"ms";
  $("#output_path").value = ST.output_path||"";
  setEnable($("#output_format"), ready && !ST.run_sfs_only);
  setEnable($("#output_path"), ready && !ST.run_sfs_only);
  setEnable($("#temp_loc"), ready);

  // engine change toggles later parts
}

function renderOptimization(){
  const ready = DV.ordered_ready;
  const showParallel = (DV.reps>1 && DV.max_workers>1);
  const parBox = $("#parallelBox");
  show(parBox, showParallel);
  const par = $("#parallel"); const parVal=$("#parVal");
  par.max = DV.max_workers||1;
  par.value = Math.min(ST.parallel||1, parseInt(par.max));
  parVal.textContent = par.value;
  setEnable(par, ready && showParallel);

  const simsBox = $("#simsBox");
  show(simsBox, ST.engine!=="msprime");
  const simsAuto = $("#sims_auto");
  simsAuto.checked = !!ST.sims_per_work_auto;
  const spw = $("#sims_per_work");
  const pval = $("#spwVal");
  const maxSpw = Math.max(1, Math.floor((DV.reps||1)/(ST.parallel||1)));
  spw.max = maxSpw; if(!spw.value) spw.value = 1;
  pval.textContent = spw.value;
  setEnable(simsAuto, ready && ST.engine!=="msprime");
  setEnable(spw, ready && ST.engine!=="msprime" && !simsAuto.checked);

  const mr = $("#max_ram_percent");
  mr.value = ST.max_ram_percent||80;
  const mrLbl = $("#ramVal");
  if(mrLbl){ mrLbl.textContent = (ST.max_ram_percent||80) + "%"; }
  setEnable(mr, ready);
  // ram cap radios
  const caps = $$('input[name="ram_cap"]');
  caps.forEach(r=>{ r.checked = (r.value === (ST.max_ram_cap||"system")); r.disabled = !ready; });

  const gmf = $("#growth_max_fold");
  gmf.value = ST.growth_max_fold || "1.05";
  show(gmf.closest(".row"), ST.engine!=="msprime");
  setEnable(gmf, ready && ST.engine!=="msprime");
}

function renderSFS(){
  const ready = DV.ordered_ready;
  const sfsOn = !!ST.sfs_on;
  const sfsChk = $("#sfs_on");
  sfsChk.checked = sfsOn; setEnable(sfsChk, ready);
  toggleHide($("#sfsPanel"), !sfsOn);
  $("#sfs_output").value = ST.sfs_output || "";
  $("#sfs_normalized").checked = !!ST.sfs_normalized;
  $("#sfs_mode").value = ST.sfs_mode || "mean";
  ["#sfs_output","#sfs_normalized","#sfs_mode"].forEach(id=> setEnable($(id), ready && sfsOn));
}

function renderSelection(){
  const ready = DV.ordered_ready;
  const supportsSweep = (ST.engine==="msms" || ST.engine==="discoal");
  const sweepChk = $("#sweep_enable");
  setEnable(sweepChk, ready && supportsSweep);
  sweepChk.checked = !!ST.sweep_enable;

  toggleHide($("#selPanel"), !(ready && supportsSweep && ST.sweep_enable));
  if(!(ready && supportsSweep && ST.sweep_enable)) { show($("#pairedCard"), false); return; }

  // sweep pop
  const sp = $("#sweep_pop");
  const options = (ST.populations||[]).map(p=>p.name);
  sp.innerHTML = `<option value="">-- select population --</option>` + options.map(n=>`<option ${ST.sweep_pop===n?'selected':''}>${n}</option>`).join("");

  // position slider: bp if we know length
  const len = ST.seq_length>0 ? ST.seq_length : (ST.chrom_length_mode && DV.contig_len ? DV.contig_len : null);
  const boxPct = $("#posPercentBox"); const boxBp = $("#posBpBox");
  if(len){
    toggleHide(boxPct, true);
    toggleHide(boxBp, false);
    const bp = $("#sweep_pos_bp");
    bp.max = len; bp.value = Math.min(ST.sweep_pos_raw||0, len);
    $("#posBpHint").textContent = `0..${Number(len).toLocaleString()} bp`;
  }else{
    toggleHide(boxPct, false);
    toggleHide(boxBp, true);
    $("#sweep_pos_percent").value = ST.sweep_pos_raw || 50;
    $("#posPctLbl").textContent = (ST.sweep_pos_raw || 50) + "%";
  }

  $("#sel_s").value = ST.sel_s || 0.1;
  $("#time_mode").value = ST.time_mode || "";
  $("#time_value").value = ST.time_mode==="Fixation Time" ? (ST.fixation_time||"") : (ST.sweep_time||"");
  $$('input[name="time_units"]').forEach(r=>{ r.checked = (r.value===(ST.time_units||"gens")); });

  // Paired neutral card visible only when sweep enabled & supported
  show($("#pairedCard"), true);
  const pnChk = $("#paired_neutral"); pnChk.checked = !!ST.paired_neutral;
  toggleHide($("#pairedPanel"), !pnChk.checked);
  $("#paired_neutral_name").value = ST.paired_neutral_name || "";
  // Render neutral engine options dynamically; show an asterisk next to
  // the engine option that matches the currently selected sweep engine.
  const neutralSel = $("#neutral_engine");
  const neutralEngines = ['discoal','ms','msms','scrm','msprime'];
  neutralSel.innerHTML = neutralEngines.map(e => {
    const mark = (e === ST.engine) ? ' *' : '';
    return `<option value="${e}" ${((ST.neutral_engine||ST.engine)===e)?'selected':''}>${e}${mark}</option>`;
  }).join('');
  // Default neutral selection to the sweep engine when none is explicitly set
  neutralSel.value = ST.neutral_engine || ST.engine;
}

function renderBuild(){
  $("#progress_flag").checked = !!ST.progress_flag;
  // Disable progress flag until species/model/chromosome/populations are configured
  setEnable($("#progress_flag"), !!DV.ordered_ready);
  const sfsOnly = $("#run_sfs_only");
  sfsOnly.checked = !!ST.run_sfs_only;
  setEnable(sfsOnly, DV.ordered_ready && ST.sfs_on);
}

function attachHandlers(){
  $("#species").addEventListener("change", e=> apiSet({species_id: e.target.value, model_id:"", chromosome:""}));
  $("#model").addEventListener("change", e=> apiSet({model_id: e.target.value, chromosome:""}));
  $("#chromosome").addEventListener("change", e=> apiSet({chromosome: e.target.value}));

  $("#engine").addEventListener("change", e=> apiSet({engine: e.target.value}));
  $("#replicates").addEventListener("change", e=> apiSet({replicates: e.target.value}));
  $("#chrom_length_mode").addEventListener("change", e=> apiSet({chrom_length_mode: e.target.checked}));
  $("#seq_length").addEventListener("change", e=> apiSet({seq_length: e.target.value}));
  $("#target_snps").addEventListener("change", e=> apiSet({target_snps: e.target.value}));
  $("#tol_auto").addEventListener("change", e=> apiSet({target_snps_tol_auto: e.target.checked}));
  $("#target_snps_tol").addEventListener("change", e=> apiSet({target_snps_tol: e.target.value}));
  $("#output_format").addEventListener("change", e=> apiSet({output_format: e.target.value}));
  $("#output_path").addEventListener("change", e=> apiSet({output_path: e.target.value}));
  $("#temp_loc").addEventListener("change", e=> apiSet({temp_loc: e.target.value}));

  $("#parallel").addEventListener("input", e=> { $("#parVal").textContent = e.target.value; });
  $("#parallel").addEventListener("change", e=> apiSet({parallel: e.target.value}));
  $("#sims_auto").addEventListener("change", e=> apiSet({sims_per_work_auto: e.target.checked}));
  $("#sims_per_work").addEventListener("input", e=> { $("#spwVal").textContent = e.target.value; });
  $("#sims_per_work").addEventListener("change", e=> apiSet({sims_per_work: e.target.value}));
  $("#max_ram_percent").addEventListener("change", e=> apiSet({max_ram_percent: e.target.value}));
  $("#max_ram_percent").addEventListener("input", e=> { const lbl = $("#ramVal"); if(lbl){ lbl.textContent = e.target.value + "%"; } });
  $$('input[name="ram_cap"]').forEach(r=> r.addEventListener("change", e=> apiSet({max_ram_cap: e.target.value})));
  $("#growth_max_fold").addEventListener("change", e=> apiSet({growth_max_fold: e.target.value}));

  $("#sfs_on").addEventListener("change", e=> apiSet({sfs_on: e.target.checked}));
  $("#sfs_output").addEventListener("change", e=> apiSet({sfs_output: e.target.value}));
  $("#sfs_normalized").addEventListener("change", e=> apiSet({sfs_normalized: e.target.checked}));
  $("#sfs_mode").addEventListener("change", e=> apiSet({sfs_mode: e.target.value}));

  $("#sweep_enable").addEventListener("change", e=> apiSet({sweep_enable: e.target.checked}));
  $("#sweep_pop").addEventListener("change", e=> apiSet({sweep_pop: e.target.value}));
  $("#sweep_pos_percent").addEventListener("input", e=> { $("#posPctLbl").textContent = e.target.value+"%"; });
  $("#sweep_pos_percent").addEventListener("change", e=> apiSet({sweep_pos_raw: e.target.value}));
  $("#sweep_pos_bp").addEventListener("change", e=> apiSet({sweep_pos_raw: e.target.value}));
  $("#sel_s").addEventListener("change", e=> apiSet({sel_s: e.target.value}));
  $("#time_mode").addEventListener("change", e=> apiSet({time_mode: e.target.value}));
  $("#time_value").addEventListener("change", e=> {
    if($("#time_mode").value==="Fixation Time"){ apiSet({fixation_time: e.target.value}); }
    else { apiSet({sweep_time: e.target.value}); }
  });
  $$('input[name="time_units"]').forEach(r=> r.addEventListener("change", e=> apiSet({time_units: e.target.value})));

  $("#paired_neutral").addEventListener("change", e=> apiSet({paired_neutral: e.target.checked}));
  $("#paired_neutral_name").addEventListener("change", e=> apiSet({paired_neutral_name: e.target.value}));
  $("#neutral_engine").addEventListener("change", e=> apiSet({neutral_engine: e.target.value}));

  $("#progress_flag").addEventListener("change", e=> apiSet({progress_flag: e.target.checked}));
  $("#run_sfs_only").addEventListener("change", e=> apiSet({run_sfs_only: e.target.checked}));

  $("#showCmd").addEventListener("click", async (e)=>{
    e.preventDefault();
    const r = await fetch("/api/show_engine", {method:"POST"});
    const js = await r.json();
    const box = $("#engineCmd");
    if(js.ok){ box.textContent = js.engine; toggleHide(box,false); }
    else {
      box.textContent = (js.errors||["Error"]).join("\\n");
      toggleHide(box,false);
    }
  });
  $("#downloadCmd").addEventListener("click", async (e)=>{
    e.preventDefault();
    const r = await fetch("/api/show_engine", {method:"POST"});
    const js = await r.json();
    if(js.ok){
      const blob = new Blob([js.engine], {type:"text/plain"});
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      // Derive a friendly filename from Primary Output base or SFS base
      const stripKnownExt = (nm)=>{
        if(!nm) return "engine_command";
        let b = nm.split("/").pop();
        const ends = [".ms.gz",".vcf.gz",".ms",".vcf",".bcf",".sfs"];
        for(const e of ends){ if(b.toLowerCase().endsWith(e)){ b = b.slice(0, -e.length); break; } }
        // final fallback remove any extension
        if(b.includes(".")) b = b.substring(0, b.lastIndexOf("."));
        return b || "engine_command";
      };
      const out = (ST.output_path||"").trim();
      const sfs = (ST.sfs_output||"").trim();
      const base = stripKnownExt(out || sfs);
      a.download = `${base}_engine_cmd.txt`;
      document.body.appendChild(a); a.click(); a.remove();
    }
  });
}

function renderAll(){
  renderSpeciesModelsContigs();
  renderPops();
  renderCore();
  renderOptimization();
  renderSFS();
  renderSelection();
  renderBuild();
  validateBuild();
}

function showErrors(errs){
  const box = document.getElementById('errors');
  if(!box) return;
  if(!errs || errs.length===0){ box.innerHTML = ''; box.style.display = 'none'; return; }
  box.style.display = '';
  if(Array.isArray(errs)){
    box.innerHTML = '<ul>' + errs.map(e=>`<li>${e}</li>`).join('') + '</ul>';
  } else if(typeof errs === 'string'){
    box.innerHTML = `<div>${errs}</div>`;
  } else {
    box.innerHTML = '';
  }
}

window.addEventListener("DOMContentLoaded", ()=>{
  renderAll();
  attachHandlers();
});

// --- Execute / Cancel / SSE wiring ---
(function(){
  const execBtn = document.getElementById('btn-execute');
  const cancelBtn = document.getElementById('btn-cancel');
  const coreBar = document.getElementById('core-progress');
  const pnBar = document.getElementById('pn-progress');
  const etaEl = document.getElementById('eta-text');
  const artifactsList = document.getElementById('artifacts-list');
  let es = null;

  function setProgress(core,pn,eta){
    if(core!=null && coreBar) coreBar.value = core;
    if(pn!=null && pnBar) pnBar.value = pn;
    if(etaEl) etaEl.textContent = eta!=null? (eta+"s") : "";
  }

  function showArtifacts(obj){
    if(!artifactsList) return;
    artifactsList.innerHTML = '';
    const keys = Object.keys(obj);
    keys.sort();
    for(const k of keys){
      const url = obj[k];
      const li = document.createElement('li');
      const a = document.createElement('a');
      // url is already /artifacts/<run>/<file>
      a.href = url;
      // Show the relative artifact path under the run directory (k)
      a.textContent = k + ': ' + url;
      li.appendChild(a);
      artifactsList.appendChild(li);
    }
  }

  if(execBtn){
    execBtn.addEventListener('click', async ()=>{
      execBtn.disabled = true; cancelBtn.disabled = false;
      // show progress/artifacts UI when a run begins
      show($('#progress-core-row'), true);
      show($('#progress-pn-row'), true);
      show($('#eta-row'), true);
      show($('#artifacts-row'), true);
      // push current state first
      try{ await fetch('/api/set', {method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify({})}); }catch(e){}
      const r = await fetch('/api/execute', {method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify({})});
  const j = await r.json();
  if(!j.ok){ showErrors(j.errors || [j.error || 'Execute failed']); execBtn.disabled=false; cancelBtn.disabled=true; return; }
  showErrors([]);
  ARTS = Object.assign({}, j.artifacts||{});
  showArtifacts(ARTS);
      es = new EventSource('/api/progress');
      es.onmessage = function(ev){ try{ const d=JSON.parse(ev.data); setProgress(d.core_pct,d.pn_pct,d.eta); }catch(e){} };
      es.addEventListener('artifact', function(ev){
        try{
          const d = JSON.parse(ev.data);
          if(d && d.artifacts){
            ARTS = Object.assign({}, ARTS, d.artifacts);
            showArtifacts(ARTS);
          }
        }catch(e){}
      });
      es.addEventListener('done', function(ev){ try{ const d=JSON.parse(ev.data); setProgress(100,100,0); }catch(e){}; execBtn.disabled=false; cancelBtn.disabled=true; es.close(); es=null; });
      es.addEventListener('error', function(ev){ console.error('SSE error', ev); });
    });
  }

  if(cancelBtn){
    cancelBtn.addEventListener('click', async ()=>{
      cancelBtn.disabled = true;
      await fetch('/api/cancel', {method:'POST'});
      if(es){ es.close(); es=null; }
      setProgress(0,0,null);
      // hide progress/artifacts UI when cancelled
      show($('#progress-core-row'), false);
      show($('#progress-pn-row'), false);
      show($('#eta-row'), false);
      show($('#artifacts-row'), false);
      // clear any errors when cancelling
      showErrors([]);
      execBtn.disabled = false; cancelBtn.disabled = true;
    });
  }
})();

// Validate build readiness and enable/disable command/run buttons
async function validateBuild(){
  const showBtn = $('#showCmd');
  const dlBtn = $('#downloadCmd');
  const execBtn = $('#btn-execute');
  try{
    const r = await fetch('/api/validate', {method:'POST'});
    const js = await r.json();
    const ok = !!js.ok;
    if(showBtn) setEnable(showBtn, ok);
    if(dlBtn) setEnable(dlBtn, ok);
    if(execBtn) setEnable(execBtn, ok && DV.ordered_ready);
    // show Streamlit-style messages when there are validation errors
    // show validation errors (priority) or derived informational messages
    if(!ok){ showErrors(js.errors || []); }
    else if(DV && Array.isArray(DV.messages) && DV.messages.length>0){ showErrors(DV.messages); }
    else { showErrors([]); }
  }catch(e){
    if(showBtn) setEnable(showBtn, false);
    if(dlBtn) setEnable(dlBtn, false);
    if(execBtn) setEnable(execBtn, false);
    showErrors(['Validation failed (network error)']);
  }
}
