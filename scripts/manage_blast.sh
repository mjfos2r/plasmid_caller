#!/bin/bash
# scripts/manage_blast.sh
set -eux

# config
BLAST_VERSION="2.16.0"
VENDOR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/vendor/blast"
BLAST_BINARIES=("blastn" "blastp" "blastx" "tblastn" "tblastx" "makeblastdb" "blastdbcmd")

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color / line termination

usage() {
    echo "Usage: $0 {check|install|path}"
    echo "  check   - Check if BLAST is installed (system or local)"
    echo "  install - Install BLAST locally"
    echo "  path    - Get path to BLAST binaries (installs if necessary)"
}

# compare two version strings
compare_versions() {
    # $1: current version $2: required version
    if [ "$(printf '%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]; then
        return 0 # 1 >= 2
    else
        return 1 # 1 < 2
    fi
}

# check if we've got a system blast executable that's compliant with the specified version.
check_system_blast() {
    echo -e "${YELLOW}Checking for system BLAST installation...${NC}"
    # check for the command and whether or not it's on path
    if command -v blastn &> /dev/null; then
        local version=$(blastn -version | head -n 1 | sed -E 's/.*([0-9]+\.[0-9]+\.[0-9]+).*/\1/')
        if [[ -n "$version" ]]; then
            echo -e "${GREEN}Found system BLAST version ${version}${NC}"

            # compare found version to required version.
            if version_compare "$version" "$BLAST_VERSION"; then
                local all_exist=true
                for binary in "${BLAST_BINARIES[@]}"; do
                    if ! command -v "$binary" &> /dev/null; then
                        all_exist=false
                        echo -e "${YELLOW}Warning: '$binary' not found in PATH${NC}"
                    fi
                done

                if $all_exist; then
                    echo -e "${GREEN}System BLAST meets requirements and will be used.${NC}"
                    return 0
                else
                    echo -e "${YELLOW}Some BLAST binaries are missing from PATH.${NC}"
                fi
            else
                echo -e "${YELLOW}System BLAST version: $version is older than required: $BLAST_VERSION.${NC}"
            fi
        else
            echo -e "${YELLOW}Could not determine system BLAST version.${NC}"
        fi
    else
        echo -e "${YELLOW}System BLAST not found in PATH${NC}"
    fi

    return 1
}

check_local_blast() {
    if [[ -d "$VENDOR_DIR/bin" ]]; then
        if [[ -f "$VENDOR_DIR/bin/blastn" ]]; then
            echo -e "${GREEN}Local BLAST installation found in $VENDOR_DIR${NC}"
            return 0
        fi
    fi
    echo -e "${YELLOW}Local BLAST installation not found${NC}"
    return 1
}

get_download_url() {
    local os=$(uname -s)
    local arch=$(uname -m)
    local url=""

    case "$os" in
        "Linux")
            if [[ "$arch" == "x86_64" ]]; then
                url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-x64-linux.tar.gz"
            elif [[ "$arch" == "aarch64" ]] || [[ "$arch" == "arm64" ]]; then
                url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-arm64-linux.tar.gz"
            else
                echo -e "${RED}Unsupported architecture: $arch${NC}"
            fi
            ;;
        "Darwin")
            url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VERSION/ncbi-blast-$BLAST_VERSION+-universal-macosx.tar.gz"
            ;;
        *)
            echo -e "${RED}Unsupported operating system: $os${NC}"
            ;;
    esac
    echo "$url"
}

install_local_blast() {
    echo -e "${YELLOW}Installing BLAST version $BLAST_VERSION locally...${NC}"
    mkdir -p "$VENDOR_DIR"

    local download_url=$(get_download_url)
    local download_md5="${download_url}.md5"
    local filename=$(basename "$download_url")
    local download_path="$VENDOR_DIR/$filename"
    local download_md5_path="$VENDOR_DIR/${filename}.md5"

    echo -e "${YELLOW}Downloading from: $download_url${NC}"
    echo -e "${YELLOW}Saving output to: $download_path${NC}"

    if command -v curl &>/dev/null; then
        curl -L "$download_url" -o "$download_path"
        curl -L "$download_md5" -o "$download_md5_path"
    elif command -v wget &>/dev/null; then
        wget "$download_url" -O "$download_path"
        wget "$download_md5" -O "$download_md5_path"
    else
        echo -e "${RED}ERROR: Need either curl or wget to proceed. ensure those are installed and try again!${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Validating md5sum...${NC}"
    MD5SUM_STATUS=$(md5sum -c "$download_md5_path" | cut -d" " -f2)
    if [[ "$MD5SUM_STATUS" != "OK" ]]; then
        echo -e "${RED}ERROR: Corrupted download. Attempting to redownload."
        exit 1
    else
        echo -e "${GREEN}Checksum is Valid!${NC}"
        local temp_dir="$VENDOR_DIR/temp"
        mkdir -p "${temp_dir}"
        echo -e "${YELLOW}Extracting tar.gz archive...${NC}"
        tar -xzvf "$download_path" -C "$temp_dir" --strip-components=1
    fi

    echo -e "${YELLOW}Moving binaries from tmp to final directory.${NC}"
    cp -r "$temp_dir"/* "$VENDOR_DIR"/
    echo -e "${YELLOW}Cleaning up temp files!${NC}"
    rm -f "$download_path"
    rm -f "$download_md5_path"
    rm -rf "$temp_dir"
    chmod +x "$VENDOR_DIR/bin"/*
    echo "NCBI BLAST+ $BLAST_VERSION" > "$VENDOR_DIR/VERSION"
    echo -e "${GREEN}BLAST $BLAST_VERSION installed successfully in $VENDOR_DIR${NC}"
}

get_blast_path() {
    if check_system_blast; then
        local blastn_path=$(which blastn)
        echo "$(dirname "$blastn_path")"
        echo "system" > "$VENDOR_DIR/BLAST_SOURCE"
        return 0
    elif check_local_blast; then
        echo "$VENDOR_DIR/bin"
        echo "local" > "$VENDOR_DIR/BLAST_SOURCE"
        return 0
    else
        install_local_blast
        echo "$VENDOR_DIR/bin"
        echo "local" > "$VENDOR_DIR/BLAST_SOURCE"
        return 0
    fi
}

# main execution

if [[ $# == 1 ]]; then
    usage
    exit 1
fi

MODE=$1
case $MODE in
    "check")
        if check_system_blast; then
            exit 0
        elif check_local_blast; then
            exit 0
        else
            exit 1
        fi
        ;;
    "install")
        install_local_blast
        ;;
    "path")
        get_blast_path
        ;;
    *)
        usage
        exit 1
        ;;
esac

exit 0