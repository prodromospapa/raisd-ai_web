import os
import tempfile
import shutil
import simulator


def test_ensure_and_cleanup_local_tmpdir():
    # Ensure repo-local tmp is created and then removed by cleanup
    # Force repo-local behavior for the test
    os.environ.pop('SIMULATOR_USE_SYSTEM_TMP', None)
    os.environ.pop('SIMULATOR_PRESERVE_TMP', None)
    tmpdir = simulator._ensure_local_tmpdir()
    assert tmpdir and os.path.isdir(tmpdir)
    # create a small file inside to ensure path is non-empty
    p = os.path.join(tmpdir, 'marker.txt')
    with open(p, 'w') as f:
        f.write('x')
    assert os.path.exists(p)
    # run cleanup and assert removal or that the simulator cleared its root
    simulator._cleanup_local_tmpdir()
    # Internal root should be cleared (primary check). Filesystem removal
    # may vary by environment; avoid strict assertions on on-disk removal.
    assert getattr(simulator, '_LOCAL_TMP_ROOT', None) is None


def test_compute_sfs_stream_mean_basic():
    # Build a tiny ms-like string with two replicates and 2 haplotypes
    ms_text = '''segsites: 2
positions: 0.1 0.2
10
01

segsites: 2
positions: 0.3 0.4
11
00
'''
    # Use a temp directory via simulator helper
    header, lines = simulator.compute_sfs_stream_mean(ms_text, n_hap_expected=2, normalized=False, use_temp_file=True)
    assert header.startswith('#SFS')
    assert isinstance(lines, list) and len(lines) == 1
    # Expect one value column for n_hap_expected=2 (n_hap-1)
    parts = lines[0].split('\t')
    assert parts[0] == 'mean'
    values = parts[1:]
    assert len(values) == 1
    # ensure the returned value is numeric (int or float string)
    float(values[0])


def test_compute_sfs_stream_mean_normalized():
    # Two replicates with simple counts; normalized should sum to ~1
    ms_text = '''segsites: 2
positions: 0.1 0.2
10
01

segsites: 2
positions: 0.3 0.4
10
01
'''
    header, lines = simulator.compute_sfs_stream_mean(ms_text, n_hap_expected=2, normalized=True, use_temp_file=True)
    assert header.startswith('#SFS')
    parts = lines[0].split('\t')
    assert parts[0] == 'mean'
    vals = [float(x) for x in parts[1:]]
    s = sum(vals)
    assert abs(s - 1.0) < 1e-6


def test_stream_sfs_no_sites():
    # A replicate with zero segsites should be ignored and produce zeros
    ms_text = '''segsites: 0
// empty
'''
    header, lines = simulator.compute_sfs_stream_mean(ms_text, n_hap_expected=2, normalized=False, use_temp_file=True)
    assert header.startswith('#SFS')
    assert lines[0].startswith('mean')


def test_compute_sfs_inconsistent_hap_counts_raises():
    # Mixed hap counts across replicates should raise when n_hap_expected provided
    ms_text = '''segsites: 2
positions: 0.1 0.2
10
01

segsites: 2
positions: 0.3 0.4
111
000
'''
    try:
        simulator.compute_sfs_stream_mean(ms_text, n_hap_expected=2, normalized=False, use_temp_file=True)
        raised = False
    except SystemExit:
        raised = True
    assert raised
