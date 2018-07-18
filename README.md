# seqbot
Scripts for sequencing automation

### `watch_flexo.py`

Currently runs on local hardawre and uploads BCL files to AWS S3

### `demuxer.py`

Very similar to the `watch_flexo` script but this one is designed to kick off a demux on our local hardware, based on downloading the sample-sheet from S3.

