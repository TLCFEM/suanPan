#!/bin/bash

CURRENT_PATH="$(dirname "$(readlink -f "$0")")"

TARGET_PATH="$HOME/.local/bin"

if [[ $# == 1 ]] && [[ $1 == "--create-link" ]]; then
  echo "This script does the following tasks:"
  echo "    1. Create a link to suanPan in ~/.local/bin"
  echo "    2. Check Sublime Text and add corresponding files to ~/.config/sublime-text/Packages/User"
  echo "    3. Add a desktop file to ~/.local/share/applications"
  echo "To uninstall, remove those files manually."
  echo ""

  # make sure the link does not exist
  if [[ -f "$TARGET_PATH/suanpan" ]]; then
    echo "$TARGET_PATH/suanpan exists, now deleting it."
    echo ""
    rm "$TARGET_PATH/suanpan"
  fi

  # make sure the path exists
  if ! [[ -d $TARGET_PATH ]]; then
    mkdir -p "$TARGET_PATH"
  fi

  # create the link
  ln -s "$CURRENT_PATH/suanPan.sh" "$TARGET_PATH/suanpan"

  echo "$TARGET_PATH/suanpan is successfully created."
  echo ""

  INSTALLED=false

  # check if sublime text is installed
  ST_DIR="$HOME/.config/sublime-text-3/Packages/User"
  if [[ -d $ST_DIR ]]; then
    INSTALLED=true
  else
    ST_DIR="$HOME/.config/sublime-text/Packages/User"
    if [[ -d $ST_DIR ]]; then
      INSTALLED=true
    fi
  fi

  # install the files
  if [[ $INSTALLED == true ]]; then
    echo "{\"cmd\":[\"suanpan\"\,\"-f\"\,\"\$file\"]\,\"selector\":\"source.supan\"\,\"file_patterns\":[\"*.supan\"\,\"*.sp\"]}" >suanPan.sublime-build
    cp suanPan.sublime* "$ST_DIR"
    echo "Sublime Text installed, configuration files are copied to default folder $ST_DIR."
    echo ""
  fi

  # make sure the path exists
  if ! [[ -d "$HOME/.local/share/applications" ]]; then
    mkdir -p "$HOME/.local/share/applications"
  fi

  # desktop file
  echo -e "[Desktop Entry]\nExec=suanpan\nVersion=2.0\nType=Application\nIcon=$CURRENT_PATH\suanPan.svg\nCategories=Science\nName=suanPan\nTerminal=true\n" >"$TARGET_PATH/../share/applications/suanPan.desktop"
  echo "$HOME/.local/share/applications/suanPan.desktop is successfully created."
else
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRENT_PATH/../lib
  "$CURRENT_PATH/suanPan" "$@"
fi
