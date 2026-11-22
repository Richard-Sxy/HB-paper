#!/usr/bin/env bash
set -e
set -o pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source ./env.sh
mkdir -p "$FOAM_USER_LIBBIN" "$FOAM_USER_APPBIN" Make/logs
# Generated dependency build tree only; the original ICSFoam src remains untouched.
python3 - <<'PY'
from pathlib import Path
import shutil,os
app=Path(os.environ['DBNS_IHB_DIR']).resolve()
original=app.parents[2]/'src'
target=Path(os.environ['ICSFOAM_BUILD_SRC']).resolve()
assert target==app/'Make/ICSFoam'
for p in original.rglob('*'):
    rel=p.relative_to(original)
    if 'lnInclude' in rel.parts or any(x.startswith('linux') for x in rel.parts): continue
    if p.is_file():
        q=target/rel;q.parent.mkdir(parents=True,exist_ok=True)
        if not q.exists() or q.read_bytes()!=p.read_bytes(): shutil.copy2(p,q)
shutil.copy2(app/'compat/blockFvMatrixUpdateMatrixInterfaces.C',
             target/'blockFvMatrix/blockFvMatrixUpdateMatrixInterfaces.C')
# Application-local library overrides (same relative paths as src/).
overlay=app/'compat/ICSFoam'
for p in overlay.rglob('*'):
    if p.is_file():
        q=target/p.relative_to(overlay)
        if not q.exists() or q.read_bytes()!=p.read_bytes(): shutil.copy2(p,q)
PY
wmake -j "${BUILD_JOBS:-2}" libso "$ICSFOAM_BUILD_SRC" 2>&1 | tee Make/logs/library.log
wmake 2>&1 | tee Make/logs/application.log
