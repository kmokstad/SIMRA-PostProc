#!/bin/bash

function clone_ifem {
  # Clone IFEM
  if ! test -d ${WORKSPACE}/deps/IFEM
  then
    pushd .
    mkdir -p $WORKSPACE/deps/IFEM
    cd $WORKSPACE/deps/IFEM
    git init .
    git remote add origin https://github.com/OPM/IFEM
    git fetch --depth 1 origin $IFEM_REVISION:branch_to_build
    test $? -eq 0 || exit 1
    git checkout branch_to_build
    popd
  fi
}

# Upstreams and revisions
declare -a upstreams

declare -A upstreamRev

# Sidestreams and revisions
declare -a sidestreams

declare -A sidestreamRev

# Downstreams and revisions
declare -a downstreams

declare -A downstreamRev

IFEM_REVISION=master
if grep -qi "ifem=" <<< $ghprbCommentBody
then
  IFEM_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ifem=([0-9]+).*/\1/g'`/merge
fi

clone_ifem

source $WORKSPACE/deps/IFEM/jenkins/build-ifem-module.sh

parseRevisions
printHeader SIMRA-PostProc

build_module_and_upstreams SIMRA-PostProc

test $? -eq 0 || exit 1
