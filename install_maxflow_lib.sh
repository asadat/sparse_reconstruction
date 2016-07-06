if [ ! -d "maxflow" ]; then
	wget http://vision.csd.uwo.ca/code/maxflow-v3.01.zip
	mkdir maxflow 
	unzip maxflow-v3.01.zip -d maxflow/.
	cp maxflow_cmake.txt maxflow/CMakeLists.txt
fi
