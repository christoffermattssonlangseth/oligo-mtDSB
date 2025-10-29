#!/usr/bin/env python
# coding: utf-8

# In[16]:


import os, re, numpy as np
from PIL import Image
import imageio.v3 as iio

def natural_key(s):
    m = re.search(r'(\d+)', os.path.basename(s))
    return int(m.group(1)) if m else s

def load_and_pad(paths, bg="white", scale=0.5):
    imgs, sizes = [], []
    for p in paths:
        im = Image.open(p).convert("RGB")
        # ↓ downscale early for smaller memory footprint
        if scale != 1.0:
            w, h = im.size
            im = im.resize((int(w * scale), int(h * scale)), Image.Resampling.LANCZOS)
        imgs.append(im)
        sizes.append(im.size)

    max_w = max(w for w, h in sizes)
    max_h = max(h for w, h in sizes)

    padded = []
    for im in imgs:
        w, h = im.size
        canvas = Image.new("RGB", (max_w, max_h), color=bg)
        canvas.paste(im, ((max_w - w)//2, (max_h - h)//2))
        padded.append(np.asarray(canvas))
    return padded

# === CONFIG ===
img_dir  = "../results/rbd_highlight_plots_0.5"
gif_path = os.path.join(img_dir, "rbd_domain_0.5_highlights_small.gif")

files = [os.path.join(img_dir, f) for f in os.listdir(img_dir) if f.endswith(".png")]
files = [f for f in files if os.path.getsize(f) > 1024]
files = sorted(files, key=lambda x: int(re.search(r'cluster_(\d+)', x).group(1)))

if not files:
    raise RuntimeError(f"No PNGs found in {img_dir}")

# ↓ make it smaller and faster
frames = load_and_pad(files, bg="white", scale=0.4)

# ↓ lower resolution, faster playback (200–400 ms per frame)
iio.imwrite(
    gif_path,
    frames,
    duration=900,   # ms per frame
    loop=0,         # infinite loop
    palettesize=64, # smaller color palette → smaller file
)

print(f"✅ Saved {len(frames)} frames to {gif_path}")


# In[ ]:




