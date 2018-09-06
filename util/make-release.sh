#!/bin/bash -ex
COMMIT="$1"
RELEASE_ID="$2"
RELEASE_MSG="$3"

if [ -z "${GITHUB_TOKEN}" ]; then
	echo "Please set GitHub Token"
  exit 1
fi

function hasCommand() {
	command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand github-release
hasCommand echo
hasCommand date
hasCommand wget

# download CI builds
wget -O "plass-static_linux.tar.gz" "https://mmseqs.com/archive/${COMMIT}/plass-static_sse41.tar.gz"
wget -O "plass-static-osx.tar.gz" "https://mmseqs.com/archive/${COMMIT}/plass-osx-static_sse41.tar.gz" 
wget -O "README" "https://raw.githubusercontent.com/soedinglab/plass/master/README.md"

# create release tag 
git tag "${RELEASE_ID}" && git push --tags

# create a formal release
github-release release \
    --user soedinglab \
    --repo plass \
    --tag "${RELEASE_ID}" \
    --name "Plass Release $RELEASE_ID" \
    --description "$RELEASE_MSG" \
    --pre-release

# upload Linux static binary
github-release upload \
    --user soedinglab \
    --repo plass \
    --tag "${RELEASE_ID}" \
    --name "Plass-Linux.tar.gz" \
    --file plass-static_linux.tar.gz

# upload OSX static binary
github-release upload \
    --user soedinglab \
    --repo plass \
    --tag "${RELEASE_ID}" \
    --name "Plass-Osx.gz" \
    --file plass-static-osx.tar.gz

github-release upload \
    --user soedinglab \
    --repo plass \
    --tag "${RELEASE_ID}" \
    --name "README" \
    --file README

