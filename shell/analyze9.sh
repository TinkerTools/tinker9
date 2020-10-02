#!/bin/bash

DIR="$(cd "$(dirname $(readlink -f "${BASH_SOURCE[0]}"))" >/dev/null 2>&1 && pwd)"
${DIR}/tinker9 analyze "$@"
