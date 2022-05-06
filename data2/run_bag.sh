gnome-terminal -- roscore
sleep 2s
gnome-terminal --title="uav" --geometry="80x25+10+10" -- rosbag play uav6.bag
gnome-terminal --title="vicon" --geometry="80x25+10+160" -- rosbag play vicon6.bag

