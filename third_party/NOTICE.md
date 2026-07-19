# Third-Party Notices

## TopoToolbox Lite

`topotoolbox_lite/` is a curated collection of unmodified source files from
[TopoToolbox v2.4](https://github.com/wschwanghart/topotoolbox/tree/2.4),
base tag commit `cdc21040d86bdebc4723bc28f689344db471410b`, released on
14 June 2022 by Wolfgang Schwanghart and contributors. HydroPol2D bundles
only the MATLAB classes and methods required by its terrain-processing
workflows. The selected files are listed in
[`topotoolbox_lite/MANIFEST.txt`](topotoolbox_lite/MANIFEST.txt).

TopoToolbox is distributed under the GNU General Public License, version 3.
Its license is retained verbatim in
[`topotoolbox_lite/LICENSE.txt`](topotoolbox_lite/LICENSE.txt), and original
source headers remain in every bundled MATLAB file.

HydroPol2D does not modify, reimplement, or claim authorship of the bundled
TopoToolbox algorithms.

Three methods retain later upstream revisions already used by HydroPol2D to
preserve current behavior and support modern MATLAB optimization syntax. Their
commits are recorded in `topotoolbox_lite/UPSTREAM.md` and `MANIFEST.txt`.
