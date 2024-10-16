#!/bin/bash -e
# Improved MMseqs setup script
# Constants
export MMSEQS_FORCE_MERGE=1
ARIA_NUM_CONN=8
PDB_AWS_SNAPSHOT="20240101"
UNIREF30DB="uniref30_2302"
DEFAULT_PDB_SERVER="rsync.wwpdb.org::ftp"
DEFAULT_PDB_PORT=33444

# Arguments with defaults
WORKDIR="${1:-$(pwd)}"
PDB_SERVER="${2:-$DEFAULT_PDB_SERVER}"
PDB_PORT="${3:-$DEFAULT_PDB_PORT}"
PDB_AWS_DOWNLOAD="${4:-}"

# Skip index creation with environment variable
MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-}

# Set up working directory
cd "${WORKDIR}" || { echo "Error: Could not change to directory ${WORKDIR}"; exit 1; }

# Logging functions
log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1"
}

fail() {
    echo "Error: $1" >&2
    exit 1
}

# Check for required commands
hasCommand() {
    command -v "$1" >/dev/null 2>&1
}

# Check for available download tools
detectDownloadStrategy() {
    local STRATEGY=""
    if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
    if hasCommand curl; then STRATEGY="$STRATEGY CURL"; fi
    if hasCommand wget; then STRATEGY="$STRATEGY WGET"; fi
    if [ -z "$STRATEGY" ]; then
        fail "No download tool found in PATH. Please install aria2c, curl, or wget."
    fi
    echo "$STRATEGY"
}

STRATEGY=$(detectDownloadStrategy)

# AWS check
if [ -n "${PDB_AWS_DOWNLOAD}" ] && ! hasCommand aws; then
    fail "AWS CLI not found. Please install AWS CLI if using PDB_AWS_DOWNLOAD."
fi

# Download function using available strategies
downloadFile() {
    local URL="$1"
    local OUTPUT="$2"

    set +e
    for TOOL in $STRATEGY; do
        case "$TOOL" in
            ARIA)
                aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -d "$(dirname "$OUTPUT")" -o "$(basename "$OUTPUT")" "$URL" && set -e && return 0
                ;;
            CURL)
                curl -L -o "$OUTPUT" "$URL" && set -e && return 0
                ;;
            WGET)
                wget -O "$OUTPUT" "$URL" && set -e && return 0
                ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

# Download and prepare UniRef30 ---> 103gb zip file
setupUniRef30() {
    if [ ! -f UNIREF30_READY ]; then
        log "Downloading and setting up UniRef30 database..."
        #downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/${UNIREF30DB}.tar.gz" "${UNIREF30DB}.tar.gz"
        tar xzvf "${UNIREF30DB}.tar.gz"
        mmseqs tsv2exprofiledb "${UNIREF30DB}" "${UNIREF30DB}_db"
        if [ -z "$MMSEQS_NO_INDEX" ]; then
            mmseqs createindex "${UNIREF30DB}_db" tmp1 --remove-tmp-files 1
        fi
        ln -sf ${UNIREF30DB}_db_mapping ${UNIREF30DB}_db.idx_mapping 2>/dev/null || true
        ln -sf ${UNIREF30DB}_db_taxonomy ${UNIREF30DB}_db.idx_taxonomy 2>/dev/null || true
        touch UNIREF30_READY
    else
        log "UniRef30 is already set up."
    fi
}

# Download and prepare ColabFold EnvDB
setupColabFoldEnvDB() {
    if [ ! -f COLABDB_READY ]; then
        log "Downloading and setting up ColabFold environment database..."
        downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz" "colabfold_envdb_202108.tar.gz"
        tar xzvf "colabfold_envdb_202108.tar.gz"
        mmseqs tsv2exprofiledb "colabfold_envdb_202108" "colabfold_envdb_202108_db"
        if [ -z "$MMSEQS_NO_INDEX" ]; then
            mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1
        fi
        touch COLABDB_READY
    else
        log "ColabFold environment database is already set up."
    fi
}

# Download and prepare PDB databases
setupPDB() {
    if [ ! -f PDB_READY ]; then
        log "Downloading and setting up PDB..."
        downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/pdb100_230517.fasta.gz" "pdb100_230517.fasta.gz"
        mmseqs createdb pdb100_230517.fasta.gz pdb100_230517
        if [ -z "$MMSEQS_NO_INDEX" ]; then
            mmseqs createindex pdb100_230517 tmp3 --remove-tmp-files 1
        fi
        touch PDB_READY
    else
        log "PDB is already set up."
    fi
}

# Download and prepare PDB100 for foldseek
setupPDB100FoldSeek() {
    if [ ! -f PDB100_READY ]; then
        log "Downloading and setting up PDB100 for FoldSeek..."
        downloadFile "https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb100_foldseek_230517.tar.gz" "pdb100_foldseek_230517.tar.gz"
        tar xzvf pdb100_foldseek_230517.tar.gz pdb100_a3m.ffdata pdb100_a3m.ffindex
        touch PDB100_READY
    else
        log "PDB100 for FoldSeek is already set up."
    fi
}

# Sync PDB mmCIF files
setupPDBmmCIF() {
    if [ ! -f PDB_MMCIF_READY ]; then
        log "Syncing PDB mmCIF files..."
        mkdir -p pdb/divided pdb/obsolete
        if [ -n "${PDB_AWS_DOWNLOAD}" ]; then
            aws s3 cp --no-sign-request --recursive "s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/divided/mmCIF/" pdb/divided/
            aws s3 cp --no-sign-request --recursive "s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/obsolete/mmCIF/" pdb/obsolete/
        fi
        rsync -rlpt -v -z --delete --port=${PDB_PORT} "${PDB_SERVER}/data/structures/divided/mmCIF/" pdb/divided
        rsync -rlpt -v -z --delete --port=${PDB_PORT} "${PDB_SERVER}/data/structures/obsolete/mmCIF/" pdb/obsolete
        touch PDB_MMCIF_READY
    else
        log "PDB mmCIF files are already synced."
    fi
}

# Main setup steps
setupUniRef30
#setupColabFoldEnvDB
#setupPDB
#setupPDB100FoldSeek
#setupPDBmmCIF

log "MMseqs setup completed."
