# ============================================================
# HydroPol2D – Video Generation (HPC → Browser-compatible MP4)
# ============================================================

# -y
#   Overwrite output file without prompting.

# -framerate 12
#   Defines input frame rate (animation speed).
#   12 fps = smooth scientific animation.
#   Lower → jumpy, Higher → smoother but larger file.

# -pattern_type glob
#   Enables wildcard input (*.png).

# -i "*.png"
#   Reads all PNG frames.
#   ⚠️ Frames must be correctly ordered by filename.
#   Preferred naming: frame_0001.png, frame_0002.png, ...

# -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"
#   Ensures width and height are divisible by 2.
#   Required for H.264 with yuv420p.
#   Crops at most 1 pixel → no visible quality loss.

# -c:v libx264
#   Uses H.264 codec (required for browsers and Docusaurus).

# -crf 18
#   Quality control (lower = better quality).
#   18 → high quality (recommended)
#   22–24 → smaller file for web

# -preset slow
#   Compression efficiency vs speed.
#   slow → smaller file, slower encoding
#   medium → faster, slightly larger file

# -pix_fmt yuv420p
#   Pixel format required for compatibility with Chrome/Safari.

# -movflags +faststart
#   Moves metadata to beginning of file.
#   Enables instant playback (important for web/docs).

# ============================================================
# Copy & Paste Command
# ============================================================


ml python/3.12.1
source /oak/stanford/groups/gorelick/Marcus/ffmpeg_pyenv/bin/activate

FFMPEG=$(python3 -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())")

PNG_DIR="/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Stanford_Big_Flood_10m/Outputs/Png_Animations"
OUT_MP4="${PNG_DIR}/stanford_lidar.mp4"

$FFMPEG -y -framerate 12 \
  -pattern_type glob -i "${PNG_DIR}/*.png" \
  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
  -c:v libx264 \
  -crf 18 \
  -preset slow \
  -pix_fmt yuv420p \
  -movflags +faststart \
  "$OUT_MP4"