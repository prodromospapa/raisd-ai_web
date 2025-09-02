#!/usr/bin/env bash
set -euo pipefail
# Check if Anaconda is installed
if ! command -v conda &> /dev/null; then
    echo "Anaconda is not installed. Please install Anaconda and try again."
    exit 1
fi

install_msms() {
    if command -v msms >/dev/null 2>&1; then
        echo "[msms] Already found in PATH ($(command -v msms)); skipping install." >&2
        return 0
    fi
    (
        DOWNLOAD_PAGE="https://www.mabs.at/publications/software-msms/downloads/" # retained for reference
        INSTALL_DIR="$HOME/software/msms"
        JAR_PATH="$INSTALL_DIR/msms.jar"
        WRAPPER_PATH="$INSTALL_DIR/msms"

        echo "[msms] Using fixed msms jar URL (no scrape):" >&2
        jar_url="https://www.mabs.at/fileadmin/user_upload/p_mabs/msms3.2rc-b163.jar"
        echo "[msms] Selected jar URL: $jar_url" >&2

        mkdir -p "$INSTALL_DIR"

        tmp_dl=$(mktemp)
        echo "[msms] Downloading jar..." >&2
        if command -v curl >/dev/null 2>&1; then
            curl -fSL "$jar_url" -o "$tmp_dl"
        elif command -v wget >/dev/null 2>&1; then
            wget -O "$tmp_dl" "$jar_url"
        else
            echo "[msms] ERROR: Need curl or wget to download msms jar." >&2
            exit 1
        fi

        if [[ -f "$JAR_PATH" ]]; then
            old_size=$(stat -c %s "$JAR_PATH" 2>/dev/null || echo 0)
            new_size=$(stat -c %s "$tmp_dl")
            if [[ "$old_size" != "$new_size" ]]; then
                ts=$(date +%Y%m%d-%H%M%S)
                mv "$JAR_PATH" "$JAR_PATH.$ts.bak"
                echo "[msms] Existing jar backed up to msms.jar.$ts.bak" >&2
            fi
        fi
        mv "$tmp_dl" "$JAR_PATH"
        chmod 644 "$JAR_PATH"

        cat > "$WRAPPER_PATH" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
JAR="$(dirname "$0")/msms.jar"
if [[ ! -f "$JAR" ]]; then
    echo "ERROR: msms.jar not found next to wrapper" >&2
    exit 1
fi
exec java -jar "$JAR" "$@"
EOF
        chmod 755 "$WRAPPER_PATH"

        if ! grep -q "software/msms" "$HOME/.bashrc" 2>/dev/null; then
            echo "export PATH=\"$INSTALL_DIR:$PATH\"" >> "$HOME/.bashrc"
            echo "[msms] Added $INSTALL_DIR to PATH in ~/.bashrc" >&2
        fi
        if [[ -f "$HOME/.zshrc" ]] && ! grep -q "software/msms" "$HOME/.zshrc" 2>/dev/null; then
            echo "export PATH=\"$INSTALL_DIR:$PATH\"" >> "$HOME/.zshrc"
            echo "[msms] Added $INSTALL_DIR to PATH in ~/.zshrc" >&2
        fi

        echo "[msms] Verifying signature (should not mention 'M.F. Sanner')..." >&2
        if "$WRAPPER_PATH" 2>&1 | head -n5 | grep -qi 'M.F. Sanner'; then
            echo "[msms] ERROR: Wrong program (surface MSMS) detected. Remove it and retry." >&2
            exit 1
        fi
        echo "[msms] Test run (expect usage/help or error about missing args)..." >&2
        set +e
        "$WRAPPER_PATH" 4 1 -t 10 -r 10 1000 >/dev/null 2>&1
        rc=$?
        set -e
        echo "[msms] Basic invocation exit code: $rc" >&2
        echo "[msms] Done." >&2
    ) || echo "[msms] WARNING: installation block failed" >&2
}

install_msms

# Ensure conda environment exists early (provides compilers & gsl)
if [ ! -d "$(conda info --base)/envs/raisd-ai" ]; then
    conda env create -f raisd-ai.yml -y
    sed -i '1669,1671c\            tempdir = tempfile.TemporaryDirectory()' "$(conda info --base)/envs/raisd-ai/lib/python3.8/site-packages/stdpopsim/slim_engine.py" || true
fi

install_discoal() {
    if command -v discoal >/dev/null 2>&1; then
        echo "[discoal] Already found in PATH ($(command -v discoal)); skipping install." >&2
        return 0
    fi
    (
        set -euo pipefail
        DISCOAL_DIR="$HOME/software/discoal"
        REPO_URL="https://github.com/pavlidisLab/discoal.git"
        BIN_PATH="$DISCOAL_DIR/discoal"
        mkdir -p "$HOME/software"
        if [[ -d "$DISCOAL_DIR/.git" ]]; then
            echo "[discoal] Updating existing repo..." >&2
            git -C "$DISCOAL_DIR" fetch --all --quiet || true
            git -C "$DISCOAL_DIR" pull --ff-only || echo "[discoal] git pull non-fast-forward; skipping." >&2
        else
            echo "[discoal] Cloning repository..." >&2
            git clone --depth 1 "$REPO_URL" "$DISCOAL_DIR"
        fi
        echo "[discoal] Building with conda environment toolchain..." >&2
        conda run -n raisd-ai make -C "$DISCOAL_DIR" clean >/dev/null 2>&1 || true
        conda run -n raisd-ai make -C "$DISCOAL_DIR"
        if [[ ! -x "$BIN_PATH" ]]; then
            echo "[discoal] ERROR: build did not produce executable at $BIN_PATH" >&2
            exit 1
        fi
        if ! grep -q "software/discoal" "$HOME/.bashrc" 2>/dev/null; then
            echo "export PATH=\"$DISCOAL_DIR:$PATH\"" >> "$HOME/.bashrc"
            echo "[discoal] Added $DISCOAL_DIR to PATH in ~/.bashrc" >&2
        fi
        if [[ -f "$HOME/.zshrc" ]] && ! grep -q "software/discoal" "$HOME/.zshrc" 2>/dev/null; then
            echo "export PATH=\"$DISCOAL_DIR:$PATH\"" >> "$HOME/.zshrc"
            echo "[discoal] Added $DISCOAL_DIR to PATH in ~/.zshrc" >&2
        fi
        echo "[discoal] Verifying basic run (expect usage/help)..." >&2
        set +e
        "$BIN_PATH" 2 10 100 -t 5 >/dev/null 2>&1
        rc=$?
        set -e
        echo "[discoal] Basic invocation exit code: $rc" >&2
        echo "[discoal] Done." >&2
    ) || echo "[discoal] WARNING: installation block failed" >&2
}

install_discoal

# download and compile RAiSD-AI and add it to the PATH (skip if executable present)
if command -v RAiSD-AI >/dev/null 2>&1; then
    echo "[RAiSD-AI] Already in PATH ($(command -v RAiSD-AI)); skipping." >&2
elif [ ! -d ~/software/RAiSD-AI ]; then
    mydir=$(pwd)
    cd ~
    mkdir -p software
    cd software
    
    mkdir RAiSD-AI
    cd ~/software/RAiSD-AI
    wget https://github.com/alachins/RAiSD-AI/archive/refs/heads/master.zip
    unzip master.zip && rm master.zip
    mv RAiSD-AI-master/* . && rm -r RAiSD-AI-master
    source ./compile-RAiSD-AI.sh
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI' ~/.bashrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI' >> ~/.bashrc
    fi
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI' ~/.zshrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI' >> ~/.zshrc
    fi
    
    source ~/.bashrc
    source ~/.zshrc
    cd $mydir
fi

if command -v RAiSD-AI-ZLIB >/dev/null 2>&1; then
    echo "[RAiSD-AI-ZLIB] Already in PATH ($(command -v RAiSD-AI-ZLIB)); skipping." >&2
elif [ ! -d ~/software/RAiSD-AI-ZLIB ]; then
    mydir=$(pwd)
    cd ~
    mkdir -p software
    cd software

    mkdir RAiSD-AI-ZLIB
    cd ~/software/RAiSD-AI-ZLIB
    wget https://github.com/alachins/RAiSD-AI/archive/refs/heads/master.zip
    unzip master.zip && rm master.zip
    mv RAiSD-AI-master/* . && rm -r RAiSD-AI-master
    source ./compile-RAiSD-AI-ZLIB.sh
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI-ZLIB' ~/.bashrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI-ZLIB' >> ~/.bashrc
    fi
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI-ZLIB' ~/.zshrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI-ZLIB' >> ~/.zshrc
    fi

    source ~/.bashrc
    source ~/.zshrc
    cd $mydir
fi

# -------------------------------------------------------------
# Make simulator.py executable and add its directory to PATH
# -------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMULATOR="$SCRIPT_DIR/simulator.py"
if [[ -f "$SIMULATOR" ]]; then
    chmod +x "$SIMULATOR"
    # Add repo directory to PATH (only once)
    if ! grep -Fq "$SCRIPT_DIR" "$HOME/.bashrc" 2>/dev/null; then
        echo "export PATH=\"$SCRIPT_DIR:$PATH\"" >> "$HOME/.bashrc"
        echo "[simulator] Added $SCRIPT_DIR to PATH in ~/.bashrc" >&2
    fi
    if [[ -f "$HOME/.zshrc" ]] && ! grep -Fq "$SCRIPT_DIR" "$HOME/.zshrc" 2>/dev/null; then
        echo "export PATH=\"$SCRIPT_DIR:$PATH\"" >> "$HOME/.zshrc"
        echo "[simulator] Added $SCRIPT_DIR to PATH in ~/.zshrc" >&2
    fi
    echo "[simulator] simulator.py is executable and directory ensured on PATH (new shells)." >&2
else
    echo "[simulator] WARNING: simulator.py not found at $SIMULATOR" >&2
fi
