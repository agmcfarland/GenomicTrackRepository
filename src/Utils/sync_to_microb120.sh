#!/bin/bash

LOCAL_DIR="/data/GenomicTrackRepository"
REMOTE_USER="agmcfarland"
REMOTE_HOST="microb120.med.upenn.edu"
REMOTE_DIR="/home/agmcfarland"

rsync -avz "$LOCAL_DIR" "$REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR"
