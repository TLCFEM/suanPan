#!/bin/bash

CPATH="$(dirname "$(readlink -f "$0")")"

TARGET="$HOME/.local/bin"

if [[ "$#" == 1 ]] && [[ "$1" == "--create-link" ]]; then
	if ! [[ -f "$TARGET/suanpan" ]]; then
		if ! [[ -d $TARGET ]]; then
			mkdir $TARGET
		fi
		ln -s $CPATH/suanPan.sh $TARGET/suanpan
		echo "$TARGET/suanpan is successfully created."
		INSTALLED=false
		STDIR="$HOME/.config/sublime-text-3/Packages/User"
		if [[ -d $STDIR ]]; then
			INSTALLED=true
		else
			STDIR="$HOME/.config/sublime-text/Packages/User"
			if [[ -d $STDIR ]]; then
				INSTALLED=true
			fi
		fi
		if [[ "$INSTALLED" = true ]]; then
			echo {\"cmd\":[\"suanpan\"\,\"-f\"\,\"\$file\"]\,\"selector\":\"source.supan\"\,\"file_patterns\":[\"*.supan\"\,\"*.sp\"]} > suanPan.sublime-build
			cp suanPan.sublime* $STDIR
			echo "Sublime Text installed, configuration files are copied to default folder $STDIR."
		fi
	else
		echo "$TARGET/suanpan exists, please manually delete it and run again."
	fi
else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CPATH
	$CPATH/suanPan $@
fi
