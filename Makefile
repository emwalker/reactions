default: setup

check: wget-exists
wget-exists: ; @which wget > /dev/null

setup: check
	wget http://amdc.in2p3.fr/nubase/nubtab12.asc -O db/nubtab12.asc

clean:
	@echo "removing cached data"
	rm -rf ~/.lenrmc
