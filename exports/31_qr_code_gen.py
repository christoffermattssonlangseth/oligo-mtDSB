#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().system(' pip install "qrcode[pil]"')


# In[3]:


# install once

import qrcode

# your link
url = "https://github.com/christoffermattssonlangseth/oligo-mtDSB"

# create qr code object
qr = qrcode.QRCode(
    version=1,  # size (1–40), or None for automatic
    error_correction=qrcode.constants.ERROR_CORRECT_L,  # L=7%, M=15%, Q=25%, H=30% error correction
    box_size=10,  # pixel size of each box
    border=2,     # border thickness
)
qr.add_data(url)
qr.make(fit=True)

# make the image
img = qr.make_image(fill_color="black", back_color="white")

# show or save
img.show()                  # opens in default image viewer
img.save("../assets/oligo-mtDSB.png")    # saves to file


# In[ ]:




