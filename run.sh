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
# removed the files inside output/results and output/
rm -rf output/results/*
rm -rf output/*.h5
./bin/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

echo -e "${YELLOW}Animating...${RESET}"
python src/animation.py
if [ $? -ne 0 ]; then
  echo -e "${RED}Animating failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Animating done!${RESET}"
