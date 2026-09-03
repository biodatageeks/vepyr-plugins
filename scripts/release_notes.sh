#!/usr/bin/env bash
# Print the repository-specific part of a release's notes as Markdown.
#
#   scripts/release_notes.sh <previous-tag-or-empty> <new-tag>
#
# Lists which plugin manifests changed since the previous release (or, for the
# first release, every manifest shipped) and shows how to pin the new tag from
# vepyr. The release workflow prepends this to GitHub's auto-generated
# "What's Changed" list of merged pull requests.
set -euo pipefail

previous="${1-}"
next="${2:?usage: $0 <previous-tag-or-empty> <new-tag>}"

if [ -n "$previous" ]; then
  echo "## Plugin manifests changed since ${previous}"
  echo
  changes="$(git diff --name-status "${previous}..HEAD" -- 'plugins/*/*.source.toml')"
  if [ -z "$changes" ]; then
    echo "_No manifest changes since ${previous}._"
  else
    while IFS=$'\t' read -r status path rest; do
      case "$status" in
        A)  label="added" ;;
        D)  label="removed" ;;
        M)  label="modified" ;;
        R*) label="renamed from \`${path}\`"; path="$rest" ;;
        *)  label="$status" ;;
      esac
      plugin="$(basename "$(dirname "$path")")"
      echo "- **${plugin}** — \`${path}\` (${label})"
    done <<< "$changes"
  fi
else
  echo "## Plugin manifests in this first release"
  echo
  shopt -s nullglob
  manifests=(plugins/*/*.source.toml)
  shopt -u nullglob
  if [ "${#manifests[@]}" -eq 0 ]; then
    echo "_No plugin manifests found._"
  else
    for path in "${manifests[@]}"; do
      plugin="$(basename "$(dirname "$path")")"
      echo "- **${plugin}** — \`${path}\`"
    done
  fi
fi

cat <<MD

## Using this release

Pin this tag as the \`version\` of the manifest repository when building a plugin cache:

\`\`\`python
import vepyr

vepyr.build_plugin_cache(
    plugin="alphamissense",      # any directory under plugins/
    version="${next}",
    source_path="AlphaMissense_hg38.tsv.gz",
    cache_dir="/data/116_GRCh38_merged",
    plugin_cache_root="/data/plugin_cache",
)
\`\`\`

See [Building a plugin cache](https://biodatageeks.org/vepyr/plugins/#building-a-plugin-cache).
MD
