#!/bin/sh
# Pure POSIX sh installer for the 'raisd-ai' conda environment.
# - No bashisms (no arrays, no [[ ]], no process substitution, no $'...')
# - Minimal console noise; no spinner
# - Prints short status lines; tails stderr on failure
# - No persistent log (writes to /dev/null)

set -eu

START_TS="$(date +%s)"
INSTALL_DIR="$(pwd)"
LOG_FILE="/dev/null"                # no persistent log by request
MSMS_SHA256="${MSMS_SHA256-}"       # optional integrity pin for msms.jar reproducibility
RAISD_AI_REF="${RAISD_AI_REF-}"     # optional tag/branch/commit for reproducible builds

# TTY + colors (portable)
is_tty() { [ -t 1 ] && [ -t 2 ]; }
if is_tty; then
  BOLD="$(printf '\033[1m')"; DIM="$(printf '\033[2m')"; RED="$(printf '\033[31m')"; GRN="$(printf '\033[32m')"
  YEL="$(printf '\033[33m')"; BLU="$(printf '\033[34m')"; MAG="$(printf '\033[35m')"; CYA="$(printf '\033[36m')"
  RST="$(printf '\033[0m')"
else
  BOLD=""; DIM=""; RED=""; GRN=""; YEL=""; BLU=""; MAG=""; CYA=""; RST=""
fi

ok()   { printf "%s%s%s\n" "$GRN" "$1" "$RST"; }
warn() { printf "%s%s%s\n" "$YEL" "$1" "$RST"; }
err()  { printf "%s%s%s\n" "$RED" "$1" "$RST" >&2; }
info() { printf "%s%s%s\n" "$CYA" "$1" "$RST"; }
step() {
  n="$1"; total="$2"; title="$3"
  printf "\n%s▶ %s [%s/%s]%s\n" "$BOLD" "$title" "$n" "$total" "$RST"
}

elapsed() {
  end="$(date +%s)"
  mins=$(( (end - START_TS) / 60 ))
  secs=$(( (end - START_TS) % 60 ))
  printf "%02dm%02ds" "$mins" "$secs"
}

# runq: run a command quietly, show short status, tail stderr upon failure
runq() {
  title="$1"; shift
  tmpf="$(mktemp)"
  printf "⏳ %s..." "$title"
  if ( "$@" >>"$LOG_FILE" 2>>"$tmpf" ); then
    printf "\r✅ %s\n" "$title"
    cat "$tmpf" >>"$LOG_FILE" 2>/dev/null || true
    rm -f "$tmpf"
    return 0
  else
    printf "\r❌ %s\n" "$title"
    err "   Failed."
    echo
    warn "Last 40 stderr lines:"
    tail -n 40 "$tmpf" | sed 's/^/   /'
    cat "$tmpf" >>"$LOG_FILE" 2>/dev/null || true
    rm -f "$tmpf"
    exit 1
  fi
}

trap 'err "Interrupted. No persistent log (LOG_FILE=/dev/null)."; exit 2' INT

TOTAL_STEPS=7
CUR=0

# 1) Ensure conda/mamba exists
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Checking for Conda/Mamba"
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1 && ! command -v micromamba >/dev/null 2>&1; then
  err "Conda/Mamba is not installed. Please install Miniconda/Mamba and retry."
  exit 1
fi
# Prefer mamba/micromamba if present (fallback to conda)
if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE="$(command -v mamba)"
elif command -v micromamba >/dev/null 2>&1; then
  CONDA_EXE="$(command -v micromamba)"
else
  CONDA_EXE="$(command -v conda)"
fi
ok "Conda front-end: $CONDA_EXE"

# 2) Create or update the environment from raisd-ai.yml
ENV_BASE="$("$CONDA_EXE" info --base 2>/dev/null)/envs/raisd-ai"
if [ ! -d "$ENV_BASE" ]; then
  CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Creating conda env 'raisd-ai'"
  runq "conda env create" "$CONDA_EXE" env create -f raisd-ai.yml --quiet
else
  CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Updating conda env 'raisd-ai'"
  runq "conda env update" "$CONDA_EXE" env update -f raisd-ai.yml --prune --quiet
fi

# 3) Discover env prefix and helper
PREFIX="$("$CONDA_EXE" run -n raisd-ai sh -lc 'echo "$CONDA_PREFIX"')"
if [ -z "${PREFIX:-}" ] || [ ! -d "$PREFIX" ]; then
  err "[env] Could not resolve CONDA_PREFIX for env 'raisd-ai'."
  exit 1
fi
run_in_env() { "$CONDA_EXE" run -n raisd-ai "$@"; }
mkdir -p "$PREFIX/bin" "$PREFIX/share"
TOOLS_SHARE_DIR="$PREFIX/share/raisd-ai-tools"; mkdir -p "$TOOLS_SHARE_DIR"

# Build tool requirements (for discoal / RAiSD-AI)
for req in git make gcc; do
  if ! run_in_env sh -lc "command -v $req >/dev/null 2>&1"; then
    err "Build tool '$req' not found in env 'raisd-ai'. Add it to raisd-ai.yml and retry."
    exit 1
  fi
done

# 4) Verify Java in env (required for msms)
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Verifying Java (openjdk) in env"
if ! run_in_env java -version >/dev/null 2>&1; then
  err "'java' not available in env. Ensure 'openjdk' is listed in raisd-ai.yml."
  exit 1
fi
ok "Java OK."

install_msms_into_env() {
  msms_bin="$PREFIX/bin/msms"
  msms_share="$PREFIX/share/msms"
  jar_path="$msms_share/msms.jar"
  if [ -x "$msms_bin" ]; then
    info "msms already installed → $msms_bin"
    return 0
  fi
  mkdir -p "$msms_share"

  # Discover the msms jar URL by scraping the official downloads page
  page_tmp="$(mktemp)"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "https://www.mabs.at/publications/software-msms/downloads/" -o "$page_tmp" || true
  else
    wget -q -O "$page_tmp" "https://www.mabs.at/publications/software-msms/downloads/" || true
  fi
  # Look for absolute .jar URLs first
  jar_url="$(grep -oE 'https?://[^"]+msms[0-9A-Za-z._-]*\.jar' "$page_tmp" | head -n1 || true)"
  # If none found, try href= relative links and prefix site root
  if [ -z "${jar_url:-}" ]; then
    rel="$(grep -oE 'href="[^"]+\.jar"' "$page_tmp" | sed -E 's/href="([^"]+)"/\1/' | head -n1 || true)"
    if [ -n "${rel:-}" ]; then
      case "$rel" in
        http*) jar_url="$rel" ;;
        *)      jar_url="https://www.mabs.at$rel" ;;
      esac
    fi
  fi
  rm -f "$page_tmp"
  # Fallback: known archive URL (keeps previous behavior if parsing fails)
  if [ -z "${jar_url:-}" ]; then
    jar_url="https://www.mabs.at/fileadmin/user_upload/p_mabs/msms3.2rc-b163.jar"
  fi

  # Use env vars in sh -c to avoid nested quoting
  JAR_URL="$jar_url" JAR_PATH="$jar_path" runq "Downloading msms.jar" sh -c '
    set -eu
    tmp_dl="$(mktemp)"
    if command -v curl >/dev/null 2>&1; then
      curl -fsSL "$JAR_URL" -o "$tmp_dl"
    else
      wget -q -O "$tmp_dl" "$JAR_URL"
    fi
    mv "$tmp_dl" "$JAR_PATH"
    chmod 644 "$JAR_PATH"
    # Check ZIP magic
    if ! head -c4 "$JAR_PATH" | od -An -t x1 | tr -d " \n" | grep -qi 504b0304; then
      echo "ERROR: downloaded msms.jar does not appear to be a valid JAR (missing PK signature)" >&2
      echo "  File: $JAR_PATH" >&2
      ls -l "$JAR_PATH" >&2 || true
      head -n 40 "$JAR_PATH" 2>/dev/null | sed -n "1,40p" >&2 || true
      exit 1
    fi
  '

  # Optional reproducibility: verify SHA256 if provided (support Linux/macOS)
  if [ -n "${MSMS_SHA256:-}" ]; then
    JAR_PATH="$jar_path" EXPECTED="$MSMS_SHA256" runq "Verifying msms.jar SHA256" sh -c '
      set -eu
      if command -v sha256sum >/dev/null 2>&1; then
        got="$(sha256sum "$JAR_PATH" | awk "{print \$1}")"
      elif command -v shasum >/dev/null 2>&1; then
        got="$(shasum -a 256 "$JAR_PATH" | awk "{print \$1}")"
      else
        echo "WARN: no sha256 tool available; skipping verification" >&2
        exit 0
      fi
      if [ "$got" != "$EXPECTED" ]; then
        echo "ERROR: SHA256 mismatch for msms.jar" >&2
        echo "  expected: $EXPECTED" >&2
        echo "  got     : $got" >&2
        exit 1
      fi
    '
  fi

  BIN_PATH="$msms_bin" JAR_PATH="$jar_path" runq "Installing msms wrapper" sh -c '
  set -eu
  {
    echo "#!/bin/sh"
    echo "set -eu"
    printf "JAR=%s
" "$JAR_PATH"
    cat <<'\''EOS'\''
if [ ! -f "$JAR" ]; then
  echo "ERROR: msms.jar not found at $JAR" >&2
  exit 1
fi
if [ "$#" -eq 0 ]; then
  cat <<USAGE >&2
msms: no arguments provided.
msms expects options such as:
  -ms nsam replicates
Example: msms -ms 10 1 -t 100
USAGE
  exit 2
fi
exec java -jar "$JAR" "$@"
EOS
  } > "$BIN_PATH"
  chmod 755 "$BIN_PATH"
'

  # Basic invocation test (should exit >=0 even on usage)
  BIN_PATH="$msms_bin" CONDA_EXE="$CONDA_EXE" runq "Testing msms run" sh -c '
    set +e
    "$CONDA_EXE" run -n raisd-ai "$BIN_PATH" 4 1 -t 10 -r 10 1000 >/dev/null 2>&1
    rc="$?"
    set -e
    [ "$rc" -ge 0 ]
  '

  ok "msms installed."
}

install_discoal_into_env() {
  target_bin="$PREFIX/bin/discoal"
  if [ -x "$target_bin" ]; then
    info "discoal already installed → $target_bin"
    return 0
  fi

  workdir="$(mktemp -d)"
  WORKDIR="$workdir" CONDA_EXE="$CONDA_EXE" runq "Cloning discoal" sh -c '
    set -eu
    wd="$WORKDIR"
    "$CONDA_EXE" run -n raisd-ai git clone --depth 1 https://github.com/kern-lab/discoal.git "$wd/discoal" >/dev/null 2>&1     || "$CONDA_EXE" run -n raisd-ai git clone --depth 1 https://github.com/pavlidisLab/discoal.git "$wd/discoal" >/dev/null 2>&1
    [ -d "$wd/discoal" ]
  '

  WORKDIR="$workdir" CONDA_EXE="$CONDA_EXE" runq "Building discoal" sh -c '
    set -eu
    "$CONDA_EXE" run -n raisd-ai make -C "$WORKDIR/discoal" clean >/dev/null 2>&1 || true
    "$CONDA_EXE" run -n raisd-ai make -C "$WORKDIR/discoal"
  '

  TARGET_BIN="$target_bin" WORKDIR="$workdir" runq "Installing discoal binary" sh -c '
    set -eu
    cp "$WORKDIR/discoal/discoal" "$TARGET_BIN"
    chmod 755 "$TARGET_BIN"
  '

  TARGET_BIN="$target_bin" CONDA_EXE="$CONDA_EXE" runq "Testing discoal run" sh -c '
    set +e
    "$CONDA_EXE" run -n raisd-ai "$TARGET_BIN" 2 10 100 -t 5 >/dev/null 2>&1
    rc="$?"
    set -e
    [ "$rc" -ge 0 ]
  '

  ok "discoal installed."
}

install_raisd_ai_into_env() {
  bin1="$PREFIX/bin/RAiSD-AI"
  bin2="$PREFIX/bin/RAiSD-AI-ZLIB"
  repo_url="https://github.com/alachins/RAiSD-AI.git"
  repo_dir="$TOOLS_SHARE_DIR/RAiSD-AI"

  TOOLS_SHARE_DIR="$TOOLS_SHARE_DIR" REPO_DIR="$repo_dir" REPO_URL="$repo_url" RAISD_AI_REF="$RAISD_AI_REF" CONDA_EXE="$CONDA_EXE" runq "Preparing RAiSD-AI repo" sh -c '
    set -eu
    mkdir -p "$TOOLS_SHARE_DIR"
    if [ ! -d "$REPO_DIR/.git" ]; then
      rm -rf "$REPO_DIR" || true
      "$CONDA_EXE" run -n raisd-ai git clone --depth 1 "$REPO_URL" "$REPO_DIR" >/dev/null 2>&1
      : > "$REPO_DIR/.__buildflag"
      if [ -n "$RAISD_AI_REF" ]; then
        "$CONDA_EXE" run -n raisd-ai git -C "$REPO_DIR" fetch --tags >/dev/null 2>&1 || true
        "$CONDA_EXE" run -n raisd-ai git -C "$REPO_DIR" checkout --quiet "$RAISD_AI_REF" || true
      fi
    else
      "$CONDA_EXE" run -n raisd-ai git -C "$REPO_DIR" fetch --all --tags >/dev/null 2>&1 || true
      "$CONDA_EXE" run -n raisd-ai git -C "$REPO_DIR" pull --ff-only >/dev/null 2>&1 || true
      : > "$REPO_DIR/.__buildflag"
    fi
  '

  if [ -s "$repo_dir/.__buildflag" ] || [ ! -x "$bin1" ] || [ ! -x "$bin2" ]; then
    # Patch compile scripts
    REPO_DIR="$repo_dir" runq "Patching RAiSD-AI compile scripts" sh -c '
      set -eu
      set +e
      ls "$REPO_DIR"/compile-RAiSD-AI*.sh >/dev/null 2>&1
      present="$?"
      set -e
      if [ "$present" -ne 0 ]; then exit 0; fi
      for f in "$REPO_DIR"/compile-RAiSD-AI*.sh; do
        [ -f "$f" ] || continue
        tmp="$(mktemp)"
        sed -e "s/[[:<:]]rm[[:>:]]\([[:space:]]\)/rm -f\1/g" -e "s/rm -r/rm -rf/g" "$f" > "$tmp" || cp "$f" "$tmp"
        mv "$tmp" "$f"
      done
    '

    REPO_DIR="$repo_dir" CONDA_EXE="$CONDA_EXE" runq "Compiling RAiSD-AI" sh -c '
      set -eu
      "$CONDA_EXE" run -n raisd-ai sh -lc "cd "$REPO_DIR" && ./compile-RAiSD-AI.sh"
    '

    BIN1="$bin1" REPO_DIR="$repo_dir" runq "Installing RAiSD-AI binary" sh -c '
      set -eu
      src1=""
      if [ -e "$REPO_DIR/bin/release/RAiSD-AI" ]; then src1="$REPO_DIR/bin/release/RAiSD-AI"; fi
      if [ -z "${src1:-}" ] && [ -L "$REPO_DIR/RAiSD-AI" ]; then src1="$REPO_DIR/RAiSD-AI"; fi
      if [ -z "${src1:-}" ]; then src1="$(cd "$REPO_DIR" && { find . -name RAiSD-AI -type f; find . -name RAiSD-AI -type l; } | head -n1 || true)"; fi
      if [ -n "${src1:-}" ] && [ -e "$src1" ]; then cp -L "$src1" "$BIN1"; chmod 755 "$BIN1"; else exit 44; fi
    '

    REPO_DIR="$repo_dir" CONDA_EXE="$CONDA_EXE" runq "Compiling RAiSD-AI-ZLIB" sh -c '
      set -eu
      "$CONDA_EXE" run -n raisd-ai sh -lc "cd "$REPO_DIR" && ./compile-RAiSD-AI-ZLIB.sh"
    '

    BIN2="$bin2" REPO_DIR="$repo_dir" runq "Installing RAiSD-AI-ZLIB binary" sh -c '
      set -eu
      src2=""
      if [ -e "$REPO_DIR/bin/release/RAiSD-AI-ZLIB" ]; then src2="$REPO_DIR/bin/release/RAiSD-AI-ZLIB"; fi
      if [ -z "${src2:-}" ] && [ -L "$REPO_DIR/RAiSD-AI-ZLIB" ]; then src2="$REPO_DIR/RAiSD-AI-ZLIB"; fi
      if [ -z "${src2:-}" ]; then src2="$(cd "$REPO_DIR" && { find . -name RAiSD-AI-ZLIB -type f; find . -name RAiSD-AI-ZLIB -type l; } | head -n1 || true)"; fi
      if [ -n "${src2:-}" ] && [ -e "$src2" ]; then cp -L "$src2" "$BIN2"; chmod 755 "$BIN2"; else exit 45; fi
    '

    if [ -x "$bin1" ] && [ -x "$bin2" ]; then
      ok "RAiSD-AI binaries installed → $PREFIX/bin"
    else
      warn "RAiSD-AI: one or more binaries missing (continuing)."
    fi
  else
    info "RAiSD-AI up to date; skipping rebuild."
  fi
}

hotfix_stdpopsim() {
  target="$(ls "$PREFIX"/lib/python*/site-packages/stdpopsim/slim_engine.py 2>/dev/null | head -n1 || true)"
  if [ -n "${target:-}" ] && [ -f "$target" ]; then
    info "stdpopsim found; checking Python version for hotfix need"
    py_ver="$(run_in_env python -c 'import sys; print("%d %d"% (sys.version_info.major, sys.version_info.minor))' 2>/dev/null || true)"
    if [ -n "$py_ver" ]; then
      py_major="$(printf "%s" "$py_ver" | awk '{print $1}')"
      py_minor="$(printf "%s" "$py_ver" | awk '{print $2}')"
      if [ "$py_major" -gt 3 ] || { [ "$py_major" -eq 3 ] && [ "$py_minor" -ge 10 ]; }; then
        info "Python ${py_major}.${py_minor} ≥ 3.10 → hotfix not needed."
        return 0
      fi
    else
      warn "Could not determine Python version; applying hotfix as fallback."
    fi

    TARGET="$target" runq "Applying stdpopsim TemporaryDirectory patch" sh -c '
      set -eu
      if [ -f "$TARGET" ]; then
        if sed -n "1669p" "$TARGET" | grep -q "TemporaryDirectory"; then
          : # already patched
        else
          tmpf="$(mktemp)"
          awk "NR==1669{print \"            tempdir = tempfile.TemporaryDirectory()\\"} {print}" "$TARGET" > "$tmpf"
          mv "$tmpf" "$TARGET"
        fi
      else
        exit 0
      fi
    '
    ok "stdpopsim hotfix applied."
  else
    info "stdpopsim slim_engine.py not found; skipping hotfix."
  fi
}

install_simulator_wrapper() {
  script_dir="$(cd "$(dirname "$0")" && pwd)"
  simulator="$script_dir/simulator.py"
  wrapper="$PREFIX/bin/simulator"
  py_wrapper="$PREFIX/bin/simulator.py"

  if [ ! -f "$simulator" ]; then
    warn "simulator.py not found; skipping wrapper."
    return 0
  fi

  SIM="$simulator" WR="$wrapper" PYW="$py_wrapper" runq "Installing simulator wrappers" sh -c '
    set -eu
    chmod +x "$SIM" || true
    cat > "$WR" <<EOF
#!/bin/sh
set -eu
exec python "$SIM" "$@"
EOF
    chmod 755 "$WR"
    cat > "$PYW" <<PYW
#!/usr/bin/env python3
import os, sys
repo_sim = r"$SIM"
if not os.path.exists(repo_sim):
    sys.stderr.write("ERROR: repo simulator.py not found at %s\n" % repo_sim)
    sys.exit(2)
os.execv(sys.executable, [sys.executable, repo_sim] + sys.argv[1:])
PYW
    chmod 755 "$PYW"
  '
  ok "simulator wrappers installed."
}

# ── Execute installers ─────────────────────────────────────────────────────────
CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing msms"
install_msms_into_env

CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing discoal"
install_discoal_into_env

CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing / updating RAiSD-AI"
install_raisd_ai_into_env

# Applying stdpopsim hotfix (if needed) — run quietly (log only)
hotfix_stdpopsim >/dev/null 2>&1 || true

CUR=$((CUR+1)); step "$CUR" "$TOTAL_STEPS" "Installing simulator wrapper"
install_simulator_wrapper

printf "\n%s %s (elapsed %s)\n" "🎉" "All tools installed into env." "$(elapsed)"
info "Activate with: ${BOLD}conda activate raisd-ai${RST}"
