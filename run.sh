#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Compiling...${RESET}"
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compilation done!${RESET}"

echo -e "${YELLOW}Running...${RESET}"
./bin/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

# check if there exists a file called "config/tmp_animation.txt"
if [ ! -f "config/tmp_animation.txt" ]; then
  echo -e "${BLUE}Animation skipped.${RESET}"
else
  type_anim=$(cat config/tmp_animation.txt)
  rm config/tmp_animation.txt # remove the tmp file
  echo -e "${YELLOW}Animating...${RESET}"
  python src/animation.py $type_anim
  if [ $? -ne 0 ]; then
    echo -e "${RED}Animating failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Animating done!${RESET}"
fi
