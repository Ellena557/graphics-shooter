#pragma once
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "freeglut.lib")
#include<vector>
#include <random>

// Include GLEW
#include "GL/glew.h"

// Include GLFW
#include "GLFW/glfw3.h"


// Include GLM
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
using namespace glm;


const double PI = 3.141592653589793238463;

struct Tree {
	std::vector<GLfloat> tree;
	std::vector<GLfloat> tree_colors;
	std::vector<GLfloat> uves;
};

struct Spruce : public Tree {
	double turn_xy;
	double turn_xz;
	double scale;
	double bias_x;
	double bias_y;
	double bias_z;

	double last_shoot = 0; // Когда в последний раз стреляла.

	glm::vec3 centre;
	double radius;
};

struct Stump : public Tree {
	double birth;
};

struct Fireboll {
	glm::vec3 centre;
	double radius;
	glm::vec3 direction;
	float speed = 0.5f;
	double birth;

	std::vector<GLfloat> boll;
	std::vector<GLfloat> boll_colors;
	std::vector<GLfloat> uves;
};

struct SpruceFireboll {
	glm::vec3 centre;
	double radius;
	glm::vec3 direction;
	float speed = 1.0f;
	double birth;

	std::vector<GLfloat> boll;
	std::vector<GLfloat> boll_colors;
	std::vector<GLfloat> uves;
};

void make_tree(Tree* res, const std::vector<GLfloat>& tree, const std::vector<GLfloat>& tree_colors,
	double turn_xy, double turn_xz, double scale, double bias_x, double bias_y, double bias_z);

class MakeSpruce {
private:
	const std::vector<GLfloat> tree = {
			 0.0f, 0.8f, 0.0f,
			-0.2f, 0.4f,-0.2f,
			-0.2f, 0.4f, 0.2f,

			 0.0f, 0.8f, 0.0f,
			-0.2f, 0.4f,-0.2f,
			 0.2f, 0.4f,-0.2f,

			 0.0f, 0.8f, 0.0f,
			 0.2f, 0.4f, 0.2f,
			-0.2f, 0.4f, 0.2f,

			 0.0f, 0.8f, 0.0f,
			 0.2f, 0.4f, 0.2f,
			 0.2f, 0.4f,-0.2f,

			 0.2f, 0.4f, 0.2f,
			-0.2f, 0.4f,-0.2f,
			-0.2f, 0.4f, 0.2f,

			 0.2f, 0.4f, 0.2f,
			-0.2f, 0.4f,-0.2f,
			 0.2f, 0.4f,-0.2f,

			 0.0f, 0.6f, 0.0f,
			-0.3f, 0.0f,-0.3f,
			-0.3f, 0.0f, 0.3f,

			 0.0f, 0.6f, 0.0f,
			-0.3f, 0.0f,-0.3f,
			 0.3f, 0.0f,-0.3f,

			 0.0f, 0.6f, 0.0f,
			 0.3f, 0.0f, 0.3f,
			-0.3f, 0.0f, 0.3f,

			 0.0f, 0.6f, 0.0f,
			 0.3f, 0.0f, 0.3f,
			 0.3f, 0.0f,-0.3f,

			 0.3f, 0.0f, 0.3f,
			-0.3f, 0.0f,-0.3f,
			-0.3f, 0.0f, 0.3f,

			 0.3f, 0.0f, 0.3f,
			-0.3f, 0.0f,-0.3f,
			 0.3f, 0.0f,-0.3f,

			 0.0f, 0.2f, 0.0f,
			-0.4f,-0.4f,-0.4f,
			-0.4f,-0.4f, 0.4f,

			 0.0f, 0.2f, 0.0f,
			-0.4f,-0.4f,-0.4f,
			 0.4f,-0.4f,-0.4f,

			 0.0f, 0.2f, 0.0f,
			 0.4f,-0.4f, 0.4f,
			-0.4f,-0.4f, 0.4f,

			 0.0f, 0.2f, 0.0f,
			 0.4f,-0.4f, 0.4f,
			 0.4f,-0.4f,-0.4f,

			 0.4f,-0.4f, 0.4f,
			-0.4f,-0.4f,-0.4f,
			-0.4f,-0.4f, 0.4f,

			 0.4f,-0.4f, 0.4f,
			-0.4f,-0.4f,-0.4f,
			 0.4f,-0.4f,-0.4f,
			 // Ñòâîë
			-0.15f,-0.8f,-0.15f,
			-0.15f,-0.4f,-0.15f,
			-0.15f,-0.4f, 0.15f,

			-0.15f,-0.4f, 0.15f,
			-0.15f,-0.8f,-0.15f,
			-0.15f,-0.8f, 0.15f,

			-0.15f,-0.8f,-0.15f,
			-0.15f,-0.4f,-0.15f,
			 0.15f,-0.4f,-0.15f,

			 0.15f,-0.4f,-0.15f,
			-0.15f,-0.8f,-0.15f,
			 0.15f,-0.8f,-0.15f,

			 0.15f,-0.8f, 0.15f,
			 0.15f,-0.4f, 0.15f,
			-0.15f,-0.4f, 0.15f,

			-0.15f,-0.4f, 0.15f,
			 0.15f,-0.8f, 0.15f,
			-0.15f,-0.8f, 0.15f,

			 0.15f,-0.8f, 0.15f,
			 0.15f,-0.4f, 0.15f,
			 0.15f,-0.4f,-0.15f,

			 0.15f,-0.4f,-0.15f,
			 0.15f,-0.8f, 0.15f,
			 0.15f,-0.8f,-0.15f,

			 0.15f,-0.8f, 0.15f,
			-0.15f,-0.8f,-0.15f,
			-0.15f,-0.8f, 0.15f,

			 0.15f,-0.8f, 0.15f,
			-0.15f,-0.8f,-0.15f,
			 0.15f,-0.8f,-0.15f };

	const std::vector<GLfloat> tree_colors = {
		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  1.0f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,
		0.0f,  0.3f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f
	};

	const std::vector<GLfloat> uves = {
		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f
	};

public:
	void add_tree(Spruce* res);
};

class MakeStump {
private:
	const std::vector<GLfloat> stump = {
		//1
	   -0.15f,-0.8f,-0.15f,
	   -0.15f,-0.5f,-0.15f,
	   -0.15f,-0.5f, 0.15f,

	   -0.15f,-0.5f, 0.15f,
	   -0.15f,-0.8f,-0.15f,
	   -0.15f,-0.8f, 0.15f,
	   //2
	   -0.15f,-0.8f,-0.15f,
	   -0.15f,-0.5f,-0.15f,
		0.15f,-0.3f,-0.15f,

		0.15f,-0.3f,-0.15f,
	   -0.15f,-0.8f,-0.15f,
		0.15f,-0.8f,-0.15f,
		//3
		0.15f,-0.8f, 0.15f,
		0.15f,-0.3f, 0.15f,
	   -0.15f,-0.5f, 0.15f,

	   -0.15f,-0.5f, 0.15f,
		0.15f,-0.8f, 0.15f,
	   -0.15f,-0.8f, 0.15f,
	   //4
		0.15f,-0.8f, 0.15f,
		0.15f,-0.3f, 0.15f,
		0.15f,-0.3f,-0.15f,

		0.15f,-0.3f,-0.15f,
		0.15f,-0.8f, 0.15f,
		0.15f,-0.8f,-0.15f,
		//0
		0.15f,-0.8f, 0.15f,
	   -0.15f,-0.8f,-0.15f,
	   -0.15f,-0.8f, 0.15f,

		0.15f,-0.8f, 0.15f,
	   -0.15f,-0.8f,-0.15f,
		0.15f,-0.8f,-0.15f,
		// Êîëüöà
		// 1
		0.15f,-0.3105f, 0.15f,
	   -0.15f,-0.5105f,-0.15f,
	   -0.15f,-0.5105f, 0.15f,

		0.15f,-0.3105f, 0.15f,
	   -0.15f,-0.5105f,-0.15f,
		0.15f,-0.3105f,-0.15f,
		// 2
		0.13f,-0.3204f, 0.13f,
	   -0.13f,-0.4904f,-0.13f,
	   -0.13f,-0.4904f, 0.13f,

		0.13f,-0.3204f, 0.13f,
	   -0.13f,-0.4904f,-0.13f,
		0.13f,-0.3204f,-0.13f,
		// 3
		0.10f,-0.3303f, 0.10f,
	   -0.10f,-0.4703f,-0.10f,
	   -0.10f,-0.4703f, 0.10f,

		0.10f,-0.3303f, 0.10f,
	   -0.10f,-0.4703f,-0.10f,
		0.10f,-0.3303f,-0.10f,
		// 4
		0.085f,-0.3402f, 0.08f,
	   -0.085f,-0.4502f,-0.08f,
	   -0.085f,-0.4502f, 0.08f,

		0.085f,-0.3402f, 0.08f,
	   -0.085f,-0.4502f,-0.08f,
		0.085f,-0.3402f,-0.08f,
		// 5
		0.05f,-0.3601f, 0.05f,
	   -0.05f,-0.4201f,-0.05f,
	   -0.05f,-0.4201f, 0.05f,

		0.05f,-0.3601f, 0.05f,
	   -0.05f,-0.4201f,-0.05f,
		0.05f,-0.3601f,-0.05f,
		// 6
		0.035f,-0.37f, 0.03f,
	   -0.035f,-0.40f,-0.03f,
	   -0.035f,-0.40f, 0.03f,

		0.035f,-0.37f, 0.03f,
	   -0.035f,-0.40f,-0.03f,
		0.035f,-0.37f,-0.03f
	};

	const std::vector<GLfloat> stump_colors = {
		//4
		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		//3
		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		//2
		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		//1
		0.2f,  0.2f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		//0
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,

		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		0.2f,  0.2f,  0.0f,
		// Êîëüöà
		// 1
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		// 2
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,

		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		// 3
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		// 4
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,

		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		// 5
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,

		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		0.6f,  0.3f,  0.0f,
		// 6
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,

		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f,
		0.8f,  0.8f,  0.0f
	};

	const std::vector<GLfloat> uves = {
		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f
	};

public:
	void add_tree(Stump* res, const Spruce& spruce);
};

const int N_PHI = 50;
const int N_PSI = 50;

const int N_PHI_big = 100;
const int N_PSI_big = 100;

class MakeFireboll {
private:
	std::vector<GLfloat> boll;
	std::vector<GLfloat> boll_colors;
	std::vector<GLfloat> uves;

public:
	MakeFireboll();
	void add_boll(glm::vec3 position, glm::vec3 direction, int scale, Fireboll* res);
	void find_position(Fireboll* res);
};

class MakeSpruceFireboll {
private:
	std::vector<GLfloat> boll;
	std::vector<GLfloat> boll_colors;

	const int scale = 40;

public:
	MakeSpruceFireboll();
	void add_boll(glm::vec3 position, glm::vec3 direction, SpruceFireboll* res);
	void find_position(SpruceFireboll* res);
};

struct Foreground {
	std::vector<GLfloat> foreground;
	std::vector<GLfloat> foreground_color;
	std::vector<GLfloat> uves;
};

class MakeForeground {
private:
	std::vector<GLfloat> foreground;
	std::vector<GLfloat> foreground_color;
	std::vector<GLfloat> uves;

public:
	void add_foreground(glm::vec3 position, glm::vec3 direction, Foreground* res);
};

struct Floor {
	std::vector<GLfloat> floor;
	std::vector<GLfloat> colors;
	std::vector<GLfloat> uves;
};

class MakeFloor {
private:
	std::vector<GLfloat> floor;
	std::vector<GLfloat> colors;
	std::vector<GLfloat> uves;

public:
	void add_floor(glm::vec3 position, int floor_range, Floor* res);
};

struct Sky {
	std::vector<GLfloat> sky;
	std::vector<GLfloat> uves;
};

class MakeSky {
private:
	std::vector<GLfloat> sky;
	std::vector<GLfloat> uves;

public:
	void add_sky(glm::vec3 position, Sky* res);
};
