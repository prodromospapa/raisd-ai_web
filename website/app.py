from flask import Flask, render_template_string, send_from_directory, url_for
import os

app = Flask(__name__, static_folder='static', template_folder='templates')


TEMPLATE = """
<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>RAiSD-AI — Launcher</title>
    <style>
      body { font-family: Inter, Arial, sans-serif; margin:0; background:#f8fafc; color:#0f172a }
      .container { max-width:900px; margin:4rem auto; padding:2rem; background:#fff; border-radius:12px; box-shadow:0 6px 24px rgba(2,6,23,0.06);} 
      header { display:flex; align-items:center; gap:16px }
      h1 { margin:0; font-size:1.5rem }
      p.lead { color:#475569 }
      .cards { display:flex; gap:12px; margin-top:1.25rem }
      .card { flex:1; padding:1rem; border-radius:10px; background:#f1f5f9 }
      .cta { display:inline-block; margin-top:1rem; padding:10px 14px; background:linear-gradient(90deg,#6366f1,#06b6d4); color:#fff; border-radius:8px; text-decoration:none; font-weight:600 }
      .muted { color:#64748b; font-size:0.95rem }
      footer { margin-top:1.5rem; color:#94a3b8; font-size:0.9rem }
      .warn { background:#fff7ed; border-left:4px solid #f59e0b; padding:10px 12px; border-radius:6px }
    </style>
  </head>
  <body>
    <div class="container">
      <header>
        <div>
          <h1>RAiSD-AI — Launcher</h1>
          <p class="lead">Quick access to the Genomic Simulator UI and artifacts.</p>
        </div>
      </header>

      <div class="cards">
        <div class="card">
          <strong>Open Streamlit UI</strong>
          <p class="muted">Launch the interactive simulator UI built with Streamlit.</p>
          <a class="cta" href="/ui" target="_blank">Open UI</a>
        </div>
        <div class="card">
          <strong>Artifacts</strong>
          <p class="muted">Browse generated artifacts (zips, SFS, VCFs) produced by runs.</p>
          <a class="cta" href="/artifacts" target="_blank">Open Artifacts</a>
        </div>
      </div>

      <section style="margin-top:1.5rem">
        <div class="warn">
          <strong>Environment check</strong>
          <div style="margin-top:6px">
            {% if missing_simulator %}
              <p class="muted">The simulator script (<code>simulator.py</code>) or RAiSD executable was not found in the workspace root. The Streamlit UI will show guidance if required files are missing.</p>
            {% else %}
              <p class="muted">Simulator script detected. You can open the UI or run simulations from the Streamlit app.</p>
            {% endif %}
          </div>
        </div>
      </section>

      <footer>
        <div>Tips: Start the Streamlit UI with <code>streamlit run ui_simulator.py</code> (recommended inside the <code>raisd-ai</code> conda env).</div>
      </footer>
    </div>
  </body>
</html>
"""


@app.route('/')
def index():
    # Quick check for simulator presence
    ws = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    sim_path = os.path.join(ws, 'simulator.py')
    missing = not os.path.exists(sim_path)
    return render_template_string(TEMPLATE, missing_simulator=missing)


@app.route('/ui')
def ui_redirect():
    # Redirect suggestion: instruct user to open Streamlit interface
    return render_template_string('<p>Please open the Streamlit UI by running <code>streamlit run ui_simulator.py</code> in the project root.</p>')


@app.route('/artifacts')
def artifacts():
    artifacts_dir = os.path.join(os.path.dirname(__file__), '..', 'artifacts')
    artifacts_dir = os.path.abspath(artifacts_dir)
    if not os.path.isdir(artifacts_dir):
        return render_template_string('<p>No artifacts directory found.</p>')
    # show a simple listing
    files = sorted(os.listdir(artifacts_dir))
    items = ''.join(f'<li><a href="/artifacts/{f}">{f}</a></li>' for f in files)
    return render_template_string(f'<h2>Artifacts</h2><ul>{items}</ul>')


@app.route('/artifacts/<path:filename>')
def artifacts_file(filename):
    artifacts_dir = os.path.join(os.path.dirname(__file__), '..', 'artifacts')
    artifacts_dir = os.path.abspath(artifacts_dir)
    if not os.path.isdir(artifacts_dir):
        return 'Not found', 404
    return send_from_directory(artifacts_dir, filename, as_attachment=True)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=True)
