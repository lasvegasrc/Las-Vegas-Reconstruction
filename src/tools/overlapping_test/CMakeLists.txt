#####################################################################################
# Set source files
#####################################################################################

set(LVR_OVERLAPPING_TEST_SRC
    Main.cpp
    Options.cpp
)

#####################################################################################
# Setup dependencies to external libraries
#####################################################################################

set(LVR_OVERLAPPING_TEST_DEPS
	lvr_static
	lvrlas_static
	lvrrply_static
	lvrslam6d_static
	${OPENGL_LIBRARIES}
	${GLUT_LIBRARIES}
	${OpenCV_LIBS}
	)

#####################################################################################
# Add executable
#####################################################################################

###### ADD YOUR CODE HERE #######
add_executable(lvr_overlapping_test ${LVR_OVERLAPPING_TEST_SRC})

target_link_libraries(lvr_overlapping_test ${LVR_OVERLAPPING_TEST_DEPS})
