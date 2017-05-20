#####################################################################################
# Set source files
#####################################################################################

set(LVR_CUDA_NORMAL_SRC
    Main.cpp
    Options.cpp
)

#####################################################################################
# Setup dependencies to external libraries
#####################################################################################

set(LVR_CUDA_NORMAL_DEPS
	lvr_static
	lvrlas_static
	lvrrply_static
	lvrslam6d_static
	lvrcuda_static
	${OPENGL_LIBRARIES}
	${GLUT_LIBRARIES}
	${OpenCV_LIBS}
	)

#####################################################################################
# Add executable
#####################################################################################

###### ADD YOUR CODE HERE #######
cuda_add_executable(lvr_cuda_normals ${LVR_CUDA_NORMAL_SRC})

target_link_libraries(lvr_cuda_normals ${LVR_CUDA_NORMAL_DEPS})
