#!/usr/bin/env bash
# Submit the ADToolbox recipe to Bioconda: fill the sha256, fork
# bioconda-recipes, add the recipe, and open a PR.
#
# Run this ONLY after the two steps that require your credentials:
#   1. Publish 1.1.15 to PyPI *from a tree that has the distutils fix*
#      (adtoolbox/core.py must use `from warnings import warn`).
#   2. Authenticate GitHub CLI:  gh auth login
#
# The script refuses to proceed if a blocker remains, so it will not push a
# recipe that would fail Bioconda CI.
#
#   bash conda-recipe/submit_bioconda.sh
set -euo pipefail

VERSION="1.1.15"
RECIPE_DIR="$(cd "$(dirname "$0")" && pwd)"
META="$RECIPE_DIR/meta.yaml"

die() { echo "ERROR: $*" >&2; exit 1; }

# --- preconditions ---------------------------------------------------------
command -v gh   >/dev/null || die "gh CLI not installed"
command -v git  >/dev/null || die "git not installed"
gh auth status >/dev/null 2>&1 || die "not logged in — run: gh auth login"

echo ">> Checking PyPI for adtoolbox $VERSION ..."
SHA=$(curl -s "https://pypi.org/pypi/adtoolbox/json" | python3 -c "
import json,sys
d=json.load(sys.stdin); v='$VERSION'
if v not in d['releases']: sys.exit('not-published')
print([f['digests']['sha256'] for f in d['releases'][v] if f['packagetype']=='sdist'][0])
") || die "$VERSION is not published on PyPI yet — publish it first"
echo "   sdist sha256: $SHA"

echo ">> Verifying the published sdist has the distutils fix ..."
TMP=$(mktemp -d)
curl -sL "https://pypi.org/packages/source/a/adtoolbox/adtoolbox-$VERSION.tar.gz" -o "$TMP/s.tgz"
tar xzf "$TMP/s.tgz" -C "$TMP"
if grep -q "from distutils" "$TMP/adtoolbox-$VERSION/adtoolbox/core.py"; then
  die "published $VERSION STILL imports distutils — republish with the fix before submitting"
fi
echo "   OK: no distutils import in published $VERSION"

# license must be set (Bioconda linter blocks LicenseRef-UNSET)
grep -q "LicenseRef-UNSET" "$META" && die "set a real license in meta.yaml (about.license + license_file) first"

# --- write the sha256 into the recipe --------------------------------------
python3 - "$META" "$SHA" <<'PY'
import re,sys
meta,sha=sys.argv[1],sys.argv[2]
s=open(meta).read()
s=re.sub(r'sha256:\s*"[0-9a-fA-F]{64}".*', f'sha256: {sha}', s, count=1)
open(meta,'w').write(s)
print(">> meta.yaml sha256 updated")
PY

# --- fork, branch, add recipe, PR ------------------------------------------
# bioconda-recipes is huge; fork first, then shallow-clone the fork (a fresh
# fork is at parity with upstream master, so no upstream fetch is needed).
GH_USER=$(gh api user --jq .login)
echo ">> Forking bioconda/bioconda-recipes to $GH_USER ..."
gh repo fork bioconda/bioconda-recipes --clone=false >/dev/null
echo ">> Waiting for the fork to be ready ..."
for i in $(seq 1 30); do
  gh repo view "$GH_USER/bioconda-recipes" >/dev/null 2>&1 && break
  sleep 2
done

WORK=$(mktemp -d)
echo ">> Shallow-cloning the fork into $WORK ..."
git clone --depth=1 "https://github.com/$GH_USER/bioconda-recipes.git" "$WORK/bioconda-recipes" --quiet
cd "$WORK/bioconda-recipes"
git checkout -b add-adtoolbox

mkdir -p recipes/adtoolbox
cp "$META" recipes/adtoolbox/meta.yaml
# The recipe references license_file: LICENSE; ship it alongside meta.yaml
# because the published sdist does not yet include one.
[ -f "$RECIPE_DIR/LICENSE" ] && cp "$RECIPE_DIR/LICENSE" recipes/adtoolbox/LICENSE
git add recipes/adtoolbox/
git commit -m "Add adtoolbox $VERSION" --quiet
git push -u origin add-adtoolbox --quiet

gh pr create --repo bioconda/bioconda-recipes \
  --title "Add adtoolbox $VERSION" \
  --body "Adds [ADToolbox](https://github.com/chan-csu/ADToolbox) $VERSION — modeling and optimization of anaerobic digestion (ADM1/e-ADM) from sequencing reads." \
  --base master

echo ">> Done. Watch the PR's Bioconda CI for build/test results."
