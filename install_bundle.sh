#!/bin/bash -i

cmd(){
	wget --no-check-certificate https://github.com/jhuapl-bio/omics_workshop/raw/main/install.sh \
		-O install.sh \
		&& bash install.sh
}


# Detect OS and open a new terminal window running the specified command
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS - using AppleScript to open Terminal and execute the command
    osascript <<EOF
	tell application "Terminal"
    do script "wget --no-check-certificate https://github.com/jhuapl-bio/omics_workshop/raw/main/install.sh -O install.sh && bash install.sh"
    activate
end tell
EOF
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux - using GNOME Terminal as an example
	# try to run a command and if it fails then open a new terminal window
    gnome-terminal -- /bin/bash -c "wget --no-check-certificate https://github.com/jhuapl-bio/omics_workshop/raw/main/install.sh -O install.sh && bash install.sh; exec bash"
else
    cmd
fi

