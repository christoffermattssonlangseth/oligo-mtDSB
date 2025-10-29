#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from pathlib import Path
import subprocess


# In[4]:


BASE_DIR = '/Users/christoffer/work/karolinska/development/oligo-mtDSB/'


# In[6]:


# --- configuration ---
NOTEBOOK_DIRS = ["notebooks-01", "notebooks-02", "notebooks-03"]
OUTPUT_DIR = BASE_DIR + "exports/"   # or "scripts"
os.mkdir(OUTPUT_DIR, exist_ok=True)

# --- find all notebooks recursively ---
notebooks = []
for nb_dir in NOTEBOOK_DIRS:
    nb_path = BASE_DIR / nb_dir
    if nb_path.exists():
        notebooks.extend(nb_path.rglob("*.ipynb"))

if not notebooks:
    print("⚠️ No notebooks found in specified directories.")
    exit(0)

print(f"📚 Found {len(notebooks)} notebooks to export...")

# --- export each notebook using nbconvert ---
for nb in notebooks:
    print(f"→ Converting: {nb.relative_to(BASE_DIR)}")
    try:
        subprocess.run(
            [
                "jupyter", "nbconvert",
                "--to", "script",
                "--output-dir", str(OUTPUT_DIR),
                str(nb)
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        print(f"⚠️ Failed: {nb.name}\n{e.stderr.decode()}")

print(f"\n✅ Done! All .py scripts written to: {OUTPUT_DIR}")


# In[ ]:




