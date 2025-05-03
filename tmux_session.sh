#!/bin/bash

SESSION=WISArD
tmux kill-session -t $SESSION

# Create a new tmux session, starting in left pane's directory
tmux new-session -d -s $SESSION -c ~/project1
tmux send-keys -t $SESSION:0 'cd /home/local1/Documents/2024_Analysis/Grouper/; clear; clear' C-m
tmux set -g mouse on
# Split the window horizontally: left and right
tmux split-window -h -t $SESSION:0 -c ~/project2
tmux send-keys -t $SESSION:0.1 'cd /mnt/hgfs/shared-2/; clear; clear' C-m

# Now split the **right pane** vertically (into top and bottom)
tmux split-window -v -t $SESSION:0.1 -c ~/project3
tmux send-keys -t $SESSION:0.2 'htop' C-m

trap "tmux kill-session -t $SESSION" EXIT
# # tmux set-option destroy-unattached on
# Attach to the session
exec tmux attach-session -t $SESSION