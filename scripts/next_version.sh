#!/usr/bin/env bash
# Print the next semantic version tag for this repository.
#
#   scripts/next_version.sh patch|minor|major   # print the bumped tag
#   scripts/next_version.sh current             # print the latest tag ("" if none)
#
# Looks up the highest existing tag of the form vMAJOR.MINOR.PATCH (any other
# tag is ignored), bumps the requested component and prints the result, e.g.
# "v0.3.0". With no such tag yet the base is v0.0.0. Used by the release
# workflow; safe to run locally to preview what the next release will be.
set -euo pipefail

bump="${1:-}"
case "$bump" in
  patch|minor|major|current) ;;
  *) echo "usage: $0 patch|minor|major|current" >&2; exit 2 ;;
esac

current="$(git tag --list 'v[0-9]*.[0-9]*.[0-9]*' \
  | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' \
  | sort -t. -k1,1V -k2,2n -k3,3n \
  | tail -n 1 || true)"
if [ "$bump" = current ]; then
  echo "$current"
  exit 0
fi
current="${current:-v0.0.0}"

IFS=. read -r major minor patch <<< "${current#v}"
case "$bump" in
  major) major=$((major + 1)); minor=0; patch=0 ;;
  minor) minor=$((minor + 1)); patch=0 ;;
  patch) patch=$((patch + 1)) ;;
esac

echo "v${major}.${minor}.${patch}"
