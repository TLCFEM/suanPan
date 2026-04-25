#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# suanPan.sh
#
# Usage:
#
# 1. Setup mode:
#      suanPan.sh --create-link
#    - Creates symbolic links in ~/.local/bin (suanpan, sp)
#    - Copies Sublime Text configuration files if Sublime Text is installed
#    - Adds a desktop entry for suanPan in ~/.local/share/applications
#
# 2. Run mode:
#      suanPan.sh [suanPan arguments...]
#    - Forwards all arguments to the suanPan executable in bin/
#    - Example: suanPan.sh -f input.supan
#
# 3. Uninstall mode:
#      suanPan.sh --uninstall
#    - Removes ~/.local/bin links (suanpan, sp)
#    - Removes Sublime Text suanPan configuration files
#    - Removes ~/.local/share/applications/suanPan.desktop
#    - Removes ~/.suanpan-history.sp
# -----------------------------------------------------------------------------

set -e

CURRENT_PATH="$(dirname "$(readlink -f "$0")")"

TARGET_PATH="$HOME/.local/bin"

if [[ $# == 1 ]] && [[ $1 == "--create-link" ]]; then
  echo "This script does the following tasks:"
  echo "    1. Create a link to suanPan in ~/.local/bin"
  echo "    2. Check Sublime Text and add corresponding files to ~/.config/sublime-text/Packages/User"
  echo "    3. Add a desktop file to ~/.local/share/applications"
  echo "To uninstall, remove those files manually."
  echo ""

  mkdir -p "$TARGET_PATH"
  rm -f "$TARGET_PATH/sp" "$TARGET_PATH/suanpan"

  ln -s "$CURRENT_PATH/suanPan.sh" "$TARGET_PATH/suanpan"
  ln -s "$CURRENT_PATH/suanPan.sh" "$TARGET_PATH/sp"

  echo "$TARGET_PATH/suanpan and $TARGET_PATH/sp are successfully created."
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
    echo "{\"cmd\":[\"suanpan\",\"-f\",\"\$file\"],\"selector\":\"source.supan\",\"file_patterns\":[\"*.supan\",\"*.sp\"]}" >share/suanPan/suanPan.sublime-build
    cp share/suanPan/suanPan.sublime* "$ST_DIR"
    echo "Sublime Text installed, configuration files are copied to default folder $ST_DIR."
    echo ""
  fi

  # make sure the path exists
  if ! [[ -d "$HOME/.local/share/applications" ]]; then
    mkdir -p "$HOME/.local/share/applications"
  fi

  # desktop file
  echo -e "[Desktop Entry]\nExec=suanpan\nVersion=2.0\nType=Application\nIcon=$CURRENT_PATH/share/icons/hicolor/scalable/apps/suanPan.svg\nCategories=Education;Science\nName=suanPan\nTerminal=true\n" >"$TARGET_PATH/../share/applications/suanPan.desktop"
  echo "$HOME/.local/share/applications/suanPan.desktop is successfully created."
elif [[ $# == 1 ]] && [[ $1 == "--uninstall" ]]; then
  echo "Uninstalling suanPan links and configuration files..."
  rm -f "$HOME/.suanpan-history.sp" "$TARGET_PATH/suanpan" "$TARGET_PATH/sp" "$HOME/.local/share/applications/suanPan.desktop"
  echo "Removed files: $TARGET_PATH/suanpan, $TARGET_PATH/sp, $HOME/.local/share/applications/suanPan.desktop"

  ST_DIRS=("$HOME/.config/sublime-text-3/Packages/User" "$HOME/.config/sublime-text/Packages/User")
  for DIR in "${ST_DIRS[@]}"; do
    if [[ -d "$DIR" ]]; then
      rm -f "$DIR"/suanPan.sublime*
      echo "Removed Sublime Text configurations in: $DIR"
    fi
  done

  echo "Uninstall complete. Please remove the archive folder manually if needed."
else
  case "$(uname -s)" in
  Darwin*)
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CURRENT_PATH/lib
    if [[ -f "$CURRENT_PATH/lib/libmimalloc.dylib" ]]; then
      export DYLD_INSERT_LIBRARIES="$CURRENT_PATH/lib/libmimalloc.dylib"
    fi
    ;;
  *)
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRENT_PATH/lib
    if [[ -f "$CURRENT_PATH/lib/libmimalloc.so" ]]; then
      export LD_PRELOAD="$CURRENT_PATH/lib/libmimalloc.so"
    fi
    ;;
  esac
  export PATH=$PATH:$CURRENT_PATH
  "$CURRENT_PATH/bin/suanPan" "$@"
fi
