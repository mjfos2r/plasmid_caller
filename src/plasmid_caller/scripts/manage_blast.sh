#!/usr/bin/env bash
# scripts/manage_blast.sh
set -euo pipefail

# config
BLAST_VERSION="2.16.0"
VENDOR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/vendor/blast"
BLAST_BINARIES=("blastn" "blastp" "blastx" "tblastn" "tblastx" "makeblastdb" "blastdbcmd")

# Colors for output
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'; NC='\033[0m'

usage() {
    echo "Usage: $0 {check|install|path}"
    echo "  check   - Check if BLAST is installed (system or local)"
    echo "  install - Install BLAST locally"
    echo "  path    - Get path to BLAST binaries (installs if necessary)"
}

# compare two version strings
compare_versions() {
    # $1: current version $2: required version
    if [[ "$(printf '%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]; then
        return 0
    else
        return 1
    fi
}

# check if we've got a system blast executable that's compliant with the specified version.
check_system_blast() {
    echo -e "${YELLOW}Checking for system BLAST installation...${NC}"
    # check for the command and whether or not it's on path
    if command -v blastn &> /dev/null; then
        local v=$(blastn -version | head -n 1 | sed -E 's/.*([0-9]+\.[0-9]+\.[0-9]+).*/\1/')
        if [[ -n "$v" ]]; then
            echo -e "${GREEN}Found system BLAST version ${v}${NC}"
            # compare found version to required version.
            if compare_versions "$v" "$BLAST_VERSION"; then
                for b in "${BLAST_BINARIES[@]}"; do
                    command -v "$b" &>/dev/null || { echo -e "${YELLOW}Missing '$b' in PATH${NC}"; }
                done
                return 0
            fi
            echo -e "${YELLOW}System BLAST version: $v is too old. (Need: $BLAST_VERSION)${NC}"
        fi
    fi
    return 1
}

check_local_blast() {
    [[ -f "$VENDOR_DIR/bin/blastn" ]] && return 0
    return 1
}

get_download_url() {
    local os=$(uname -s) arch=$(uname -m)

    case "$os" in
        "Linux")
            [[ "$arch" == "x86_64" ]] \
                && echo "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-x64-linux.tar.gz" \
                || echo "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-arm64-linux.tar.gz"
            ;;
        "Darwin")
            echo "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-universal-macosx.tar.gz"
            ;;
        *)
            echo ""
            ;;
    esac
}

verify_checksum() {
    echo -e "${YELLOW}Validating md5 checksum${NC}"
    local expected=$(cut -d' ' -f1 "$1") file=$2 actual=""
    if command -v md5sum &>/dev/null; then actual=$(md5sum "$file" | awk '{print $1}')
    elif command -v md5 &>/dev/null; then actual=$(md5 -q "$file")
    else echo -e "${YELLOW}No md5 utility, skipping check. You should fix this!${NC}"; return 0;
    fi

    [[ $expected == "$actual" ]] \
        && echo -e "${GREEN}Checksum OK!${NC}" \
        || { echo -e "${RED}Checksum invalid!${NC}"; exit 1; }
}

install_local_blast() {
    echo -e "${YELLOW}Installing BLAST version $BLAST_VERSION locally to $VENDOR_DIR...${NC}"
    mkdir -p "$VENDOR_DIR"

    local url=$(get_download_url) file="$VENDOR_DIR/$(basename "$url")" md5_url="${download_url}.md5" md5_file="$VENDOR_DIR/$(basename "$md5_url")"

    [[ -z $url ]] && { echo -e "${RED}Unsupported platform ${NC}"; exit 1; }

    if command -v curl &>/dev/null; then
        curl -L "$url" -o "$file"
        curl -L "$md5_url" -o "$md5_file"
    elif command -v wget &>/dev/null; then
        wget "$url" -O "$file"
        wget "$md5_url" -O "$md5_file"
    else
        echo -e "${RED}ERROR: Need either curl or wget to proceed. ensure those are installed and try again!${NC}"
        exit 1
    fi

    verify_checksum "$md5_file" "$file"

    local temp_dir="$VENDOR_DIR/temp"
    mkdir -p "${temp_dir}"

    tar -xzvf "$file" -C "$temp_dir" --strip-components=1 #--no-same-owner
    cp -r "$temp_dir"/* "$VENDOR_DIR"/
    chmod +x "$VENDOR_DIR/bin"/*
    rm -rf "$temp_dir" "$file" "$md5_file"
    echo "NCBI BLAST+ $BLAST_VERSION" > "$VENDOR_DIR/BLAST_VERSION"
    echo -e "${GREEN}BLAST $BLAST_VERSION installed successfully in $VENDOR_DIR${NC}"
}

record_source() {
    mkdir -p "$VENDOR_DIR"
    echo "$1" >"$VENDOR_DIR/BLAST_SOURCE"
}

get_blast_path() {
    if check_system_blast; then record_source system; echo "$(dirname "$(command -v blastn)")"
    elif check_local_blast; then record_source local; echo "$VENDOR_DIR/bin"
    else install_local_blast; record_source local; echo "$VENDOR_DIR/bin"
    fi
}

# main execution

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

MODE=$1
case $MODE in
    check)   check_system_blast || check_local_blast ;;
    install) install_local_blast ;;
    path)  get_blast_path ;;
    *)       usage exit 1 ;;
esac

exit 0