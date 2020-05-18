// Include standard headers
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "freeglut.lib")
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <vector>

// Include GLEW
#include "GL/glew.h"

// Include GLFW
#include "GLFW/glfw3.h"
GLFWwindow* window;

// Include GLM
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
using namespace glm;

#include "common/shader.hpp"
#include "common/controls.hpp"
#include "common/text2D.hpp" 
#include "common/texture.hpp" 


#include "models.hpp"  
#include "blurplayer.hpp"  

#include "GL/glut.h" 

//for me only
#include "common/controls.cpp"
#include "common/shader.cpp"
#include "common/texture.cpp" 
#include "common/text2D.cpp" 


void drawFloor(GLuint floor_texture_id, GLuint floor_texture, GLuint vertexbuffer, GLuint uvbuffer, std::vector<GLfloat> floor_buffer, std::vector<GLfloat> floor_uv) {
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, floor_buffer.size() * 4, floor_buffer.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glBufferData(GL_ARRAY_BUFFER, floor_uv.size() * 4, floor_uv.data(), GL_STATIC_DRAW);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, floor_texture);
	glUniform1i(floor_texture_id, 0);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
	);

	// 2nd attribute buffer : UVs
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glVertexAttribPointer(
		1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		2,                                // size : U+V => 2
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		(void*)0                          // array buffer offset
	);

	glDrawArrays(GL_TRIANGLES, 0, floor_buffer.size()/3);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}


void drawSky(GLuint sky_texture_id, GLuint sky_texture, GLuint vertexbuffer, GLuint uvbuffer, std::vector<GLfloat> sky_buffer, std::vector<GLfloat> sky_uv) {
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sky_buffer.size() * 4, sky_buffer.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glBufferData(GL_ARRAY_BUFFER, sky_uv.size() * 4, sky_uv.data(), GL_STATIC_DRAW);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, sky_texture);
	glUniform1i(sky_texture_id, 0);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
	);

	// 2nd attribute buffer : UVs
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glVertexAttribPointer(
		1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		2,                                // size : U+V => 2
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		(void*)0                          // array buffer offset
	);

	glDrawArrays(GL_TRIANGLES, 0, sky_buffer.size() / 3);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}


int main(void)
{
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(1024, 768, "Spruce shooter", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID1 = LoadShaders("textureVertex.vertexshader", "textureFragment.fragmentshader");

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID1, "MVP");

	GLuint vertexbuffer;
	GLuint colorbuffer;
	GLuint uvbuffer;
	GLuint floorvertexbuffer;		// to draw floor
	GLuint flooruvbuffer;
	GLuint skyvertexbuffer;			// to draw sky
	GLuint skyuvbuffer;

	glGenBuffers(1, &vertexbuffer);
	glGenBuffers(1, &colorbuffer);
	glGenBuffers(1, &uvbuffer);
	glGenBuffers(1, &floorvertexbuffer);
	glGenBuffers(1, &flooruvbuffer);
	glGenBuffers(1, &skyvertexbuffer);
	glGenBuffers(1, &skyuvbuffer);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 10000);

	// степень детализации
	int cur_detalization = 50;

	MakeSpruce tree_maker;
	MakeStump stump_maker;
	MakeSpruceFireboll tree_boll_maker{ cur_detalization, cur_detalization };
	MakeForeground foreground_maker;
	MakeFloor floor_maker;
	MakeSky sky_maker;

	int cur_fireball_maker_num = 5;

	MakeFireboll fireboll_maker1{ 5 };
	MakeFireboll fireboll_maker2{ 10 };
	MakeFireboll fireboll_maker3{ 20 };
	MakeFireboll fireboll_maker4{ 30 };
	MakeFireboll fireboll_maker5{ 40 };
	MakeFireboll fireboll_maker6{ 50 };
	MakeFireboll fireboll_maker7{ 60 };
	MakeFireboll fireboll_maker8{ 70 };
	MakeFireboll fireboll_maker9{ 80 };
	MakeFireboll fireboll_maker10{ 90 };

	std::vector<MakeFireboll> all_fireball_makers = {
		fireboll_maker1,
		fireboll_maker2,
		fireboll_maker3,
		fireboll_maker4,
		fireboll_maker5,
		fireboll_maker6,
		fireboll_maker7,
		fireboll_maker8,
		fireboll_maker9,
		fireboll_maker10
	};

	MakeFireboll fireboll_maker = all_fireball_makers[cur_fireball_maker_num];

	std::vector<Spruce> trees;
	std::vector<Stump> stumps;
	std::vector<Fireboll> firebolls;
	std::vector<SpruceFireboll> tree_bolls;
	std::vector<Foreground> foreground;
	std::vector<Floor> floor;
	std::vector<Sky> sky;
	std::vector<ExplodedBall> exploded_balls;

	std::vector<GLfloat> fl_vertex_buffer_data = {};
	std::vector<GLfloat> fl_uv_buffer_data = {};

	std::vector<GLfloat> sky_vertex_buffer_data = {};
	std::vector<GLfloat> sky_uv_buffer_data = {};

	// Возможные цвета фона
	std::vector<glm::vec4> backgrounds = { {0.7f, 1.0f, 0.7f, 0.0f},
										   {0.7f, 0.7f, 1.0f, 0.0f},
										   {1.0f, 0.7f, 0.7f, 0.0f},
										   {0.0f, 0.0f, 0.0f, 0.0f},
										   {1.0f, 1.0f, 1.0f, 0.0f} };

	int boll_scale = 20; // Радиус фаербола 1/scale
	double interval = dis(gen) / 10000.0 * 5 + 3; // Задержка до появления первого дерева
	double start = glfwGetTime();
	int old_mause_state = GLFW_RELEASE;
	int old_plus_btn_state = GLFW_RELEASE;
	int old_minus_btn_state = GLFW_RELEASE;
	int old_ball_size_plus_state = GLFW_RELEASE;
	int old_ball_size_minus_state = GLFW_RELEASE;
	int old_defence_state = GLFW_RELEASE;
	int old_background_state = GLFW_RELEASE;
	int old_texture_state = GLFW_RELEASE;
	int old_floor_texture_state = GLFW_RELEASE;
	int old_sky_texture_state = GLFW_RELEASE;
	int old_plus_detal_state = GLFW_RELEASE;
	int old_minus_detal_state = GLFW_RELEASE;
	int backgorund_id = 0;
	int player_score = 0;
	int prev_floor_build_time = -10;	// когда в прошлый раз добавляли пол
	int prev_sky_build_time = -1000;	// когда в прошлый раз добавляли небо
	float current_scene_speed = 0.5f;
	float current_explosion_speed = 0.5f;
	float current_spruce_speed = 1.0f;

	glClearColor(backgrounds[backgorund_id][0], backgrounds[backgorund_id][1], backgrounds[backgorund_id][2], backgrounds[backgorund_id][3]);
	bool defence = false; // Защищаются ли ёлки.
	double effect_start = -10; // Время начала эффекта для игрока после попадания ёлки.

	//load Texture
	GLuint Texture = loadBMP_custom("fireearth.bmp");
	GLuint TextureID = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture2 = loadBMP_custom("fireearth2.bmp");
	GLuint TextureID2 = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture3 = loadBMP_custom("fireearth3.bmp");
	GLuint TextureID3 = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture4 = loadBMP_custom("fireearth4.bmp");
	GLuint TextureID4 = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture5 = loadBMP_custom("fireearth5.bmp");
	GLuint TextureID5 = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture6 = loadBMP_custom("fireearth6.bmp");
	GLuint TextureID6 = glGetUniformLocation(programID1, "myTextureSampler");

	// floor textures
	GLuint Texture_floor = loadBMP_custom("sand.bmp");
	GLuint TextureID_floor = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture_floor_1 = loadBMP_custom("lava.bmp");
	GLuint TextureID_floor_1 = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture_floor_2 = loadBMP_custom("grass.bmp");
	GLuint TextureID_floor_2 = glGetUniformLocation(programID1, "myTextureSampler");

	//sky textures
	GLuint Texture_sky = loadBMP_custom("sky2.bmp");
	GLuint TextureID_sky = glGetUniformLocation(programID1, "myTextureSampler");

	GLuint Texture_sky_2 = loadBMP_custom("sky1.bmp");
	GLuint Texture_sky_3 = loadBMP_custom("sky4.bmp");

	// текущая текстура для снаряда
	GLuint current_texture_num = 0;

	// текущая текстура для пола
	GLuint current_floor_texture_num = 1;

	// текущая текстура для неба
	GLuint current_sky_texture_num = 0;

	std::vector<GLuint> textures = {
		Texture,
		Texture2,
		Texture3,
		Texture4,
		Texture5,
		Texture6,
	};

	std::vector<GLuint> texture_ids = {
		TextureID,
		TextureID2,
		TextureID3,
		TextureID4,
		TextureID5,
		TextureID6
	};

	std::vector<GLuint> floor_textures = {
		Texture_floor,
		Texture_floor_1,
		Texture_floor_2
	};

	std::vector<GLuint> sky_textures = {
		Texture_sky,
		Texture_sky_2,
		Texture_sky_3
	};

	glUseProgram(programID1);

	// Initialize our little text library with the Holstein font
	initText2D("Holstein.DDS");

	int alive = 1;
	
	do {
		// Включение/выключение защиты ёлок по нажатию "D"
		int new_defence_state = glfwGetKey(window, GLFW_KEY_D);
		if (new_defence_state == GLFW_RELEASE && old_defence_state == GLFW_PRESS) {
			defence = !defence;
		}
		old_defence_state = new_defence_state;

		// Изменение фона по нажатию "B"
		int new_background_state = glfwGetKey(window, GLFW_KEY_B);
		if (new_background_state == GLFW_RELEASE && old_background_state == GLFW_PRESS) {
			backgorund_id = (backgorund_id + 1) % backgrounds.size();
			glClearColor(backgrounds[backgorund_id][0], backgrounds[backgorund_id][1], backgrounds[backgorund_id][2], backgrounds[backgorund_id][3]);
		}
		old_background_state = new_background_state;

		// Изменение текстуры снаряда по нажатию "T"
		int new_texture_state = glfwGetKey(window, GLFW_KEY_T);
		if (new_texture_state == GLFW_RELEASE && old_texture_state == GLFW_PRESS) {
			current_texture_num += 1;
			current_texture_num %= textures.size();
		}
		old_texture_state = new_texture_state;

		// Изменение текстуры пола по нажатию "F"
		int new_floor_texture_state = glfwGetKey(window, GLFW_KEY_F);
		if (new_floor_texture_state == GLFW_RELEASE && old_floor_texture_state == GLFW_PRESS) {
			current_floor_texture_num += 1;
			current_floor_texture_num %= floor_textures.size();
		}
		old_floor_texture_state = new_floor_texture_state;

		// Изменение текстуры неба по нажатию "S"
		int new_sky_texture_state = glfwGetKey(window, GLFW_KEY_S);
		if (new_sky_texture_state == GLFW_RELEASE && old_sky_texture_state == GLFW_PRESS) {
			current_sky_texture_num += 1;
			current_sky_texture_num %= sky_textures.size();
		}
		old_sky_texture_state = new_sky_texture_state;

		// ускорение сцены
		int new_plus_btn_state = glfwGetKey(window, GLFW_KEY_P);
		if (new_plus_btn_state == GLFW_RELEASE && old_plus_btn_state == GLFW_PRESS) {
			for (Fireboll fireboll : firebolls) {
				if (current_scene_speed < 3.0f) {
					current_scene_speed += 0.05f;
				}
				fireboll.speed = current_scene_speed;
			}

			for (SpruceFireboll fireboll : tree_bolls) {
				if (current_spruce_speed < 3.0f) {
					current_spruce_speed += 0.2f;
				}
				fireboll.speed = current_spruce_speed;
			}

			if (current_explosion_speed < 3.0f) {
				current_explosion_speed += 0.05f;
			}
		}
		old_plus_btn_state = new_plus_btn_state;

		// замедление сцены
		int new_minus_btn_state = glfwGetKey(window, GLFW_KEY_M);
		if (new_minus_btn_state == GLFW_RELEASE && old_minus_btn_state == GLFW_PRESS) {
			for (Fireboll fireboll : firebolls) {
				if (current_scene_speed > 0.1f) {
					current_scene_speed -= 0.05f;
				}
				fireboll.speed = current_scene_speed;
			}
			for (SpruceFireboll fireboll : tree_bolls) {
				if (current_spruce_speed > 0.2f) {
					current_spruce_speed -= 0.2f;
				}
				fireboll.speed = current_spruce_speed;
			}
			if (current_explosion_speed > 0.1f) {
				current_explosion_speed -= 0.05f;
			}
		}
		old_minus_btn_state = new_minus_btn_state;

		// уменьшение детализации
		int new_minus_detal_state = glfwGetKey(window, GLFW_KEY_MINUS);
		if (new_minus_detal_state == GLFW_RELEASE && old_minus_detal_state == GLFW_PRESS) {
			if (cur_fireball_maker_num > 0) {
				cur_fireball_maker_num -= 1;

				if (cur_fireball_maker_num == 0) {
					cur_detalization = 5;
				}
				else {
					cur_detalization -= 10;
				}
			}
			fireboll_maker = all_fireball_makers[cur_fireball_maker_num];
		}
		old_minus_detal_state = new_minus_detal_state;

		// увеличение детализации
		int new_plus_detal_state = glfwGetKey(window, GLFW_KEY_EQUAL);
		if (new_plus_detal_state == GLFW_RELEASE && old_plus_detal_state == GLFW_PRESS) {
			if (cur_fireball_maker_num < all_fireball_makers.size() - 1) {
				cur_fireball_maker_num += 1;
			}
			if (cur_fireball_maker_num == 1) {
				cur_detalization = 10;
			}
			else {
				cur_detalization += 10;
			}
			fireboll_maker = all_fireball_makers[cur_fireball_maker_num];
		}
		old_plus_detal_state = new_plus_detal_state;

		// Изменение размера фаерболов по нажатию "+ или -"
		/*int new_ball_size_state = glfwGetKey(window, GLFW_KEY_M);
		if (new_ball_size_state == GLFW_RELEASE && old_ball_size_minus_state == GLFW_PRESS) {
			boll_scale = std::min(boll_scale + 5, 100);
		}
		old_ball_size_minus_state = new_ball_size_state;
		new_ball_size_state = glfwGetKey(window, GLFW_KEY_P);
		if (new_ball_size_state == GLFW_RELEASE && old_ball_size_plus_state == GLFW_PRESS) {
			boll_scale = std::max(boll_scale - 5, 1);
		}
		old_ball_size_plus_state = new_ball_size_state;
		*/

		// Удаление деревьев, в которые попали файерболы
		for (int i = 0; i < firebolls.size(); ++i) {
			for (int j = 0; j < trees.size(); ++j) {
				// Учитываем, что фаербол двигается
				double delta_time = glfwGetTime() - firebolls[i].birth;
				glm::vec3 boll_centre = firebolls[i].centre + firebolls[i].direction * (float)delta_time * firebolls[i].speed;

				double distance = std::sqrt(std::pow(boll_centre.x - trees[j].centre.x, 2) +
					std::pow(boll_centre.y - trees[j].centre.y, 2) +
					std::pow(boll_centre.z - trees[j].centre.z, 2));
				if (distance < firebolls[i].radius + trees[j].radius) {
					// Создание пня
					stumps.resize(stumps.size() + 1);
					stump_maker.add_tree(&stumps[stumps.size() - 1], trees[j]);
					// Удаление дерева
					std::swap(trees[j], trees[trees.size() - 1]);
					trees.pop_back();
					// Удаление файербола
					firebolls[i].isExploded = 1;
					//std::swap(firebolls[i], firebolls[firebolls.size() - 1]);
					//firebolls.pop_back();
					player_score += 1;
					break;
				}
			}
		}

		// Удаление файербола, создание взрыва
		for (int i = 0; i < firebolls.size(); ++i) {
			if (firebolls[i].isExploded == 1) {
				firebolls[i].explosion_time = glfwGetTime();
				exploded_balls.resize(exploded_balls.size() + 1);
				fireboll_maker.add_exploded_ball(&exploded_balls[exploded_balls.size() - 1], &firebolls[i], firebolls[i].explosion_time, current_explosion_speed, cur_detalization);
				std::swap(firebolls[i], firebolls[firebolls.size() - 1]);
				firebolls.pop_back();
			}
		}

		// Удаление старых пней (срок жизни пня - 2 минуты)
		for (int i = 0; i < stumps.size(); ++i) {
			double cur_time = glfwGetTime();
			if (cur_time - stumps[i].birth > 120) {
				std::swap(stumps[i], stumps[stumps.size() - 1]);
				stumps.pop_back();
			}
		}
		
		// Удаление старых взрывов
		for (int i = 0; i < exploded_balls.size(); ++i) {
			double cur_time = glfwGetTime();
			if (cur_time - exploded_balls[i].birth > 5) {
				std::swap(exploded_balls[i], exploded_balls[exploded_balls.size() - 1]);
				exploded_balls.pop_back();
			}
		}

		// Добавление новых файерболов
		int new_mause_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
		if (new_mause_state == GLFW_RELEASE && old_mause_state == GLFW_PRESS && alive == 1) {
			firebolls.resize(firebolls.size() + 1);
			fireboll_maker.add_boll(getPosition(), getDirection(), boll_scale, &firebolls[firebolls.size() - 1]);
		}
		old_mause_state = new_mause_state;

		// Добавление новых ёлочных снарядов
		if (defence && alive == 1) {
			for (int i = 0; i < trees.size(); ++i) {
				glm::vec3 player_position = getPosition();
				double distance = std::sqrt(std::pow(player_position.x - trees[i].centre.x, 2) +
					std::pow(player_position.y - trees[i].centre.y, 2) +
					std::pow(player_position.z - trees[i].centre.z, 2));
				double shoot_time = glfwGetTime();
				if (distance < trees[i].radius * 100 && shoot_time - trees[i].last_shoot >= 5) {
					tree_bolls.resize(tree_bolls.size() + 1);
					tree_boll_maker.add_boll(trees[i].centre, player_position - trees[i].centre, &tree_bolls[tree_bolls.size() - 1]);
					trees[i].last_shoot = glfwGetTime();
				}
			}
		}

		// Удаление попавших ёлочных снарядов
		for (int i = 0; i < tree_bolls.size(); ++i) {
			// Учитываем, что фаербол двигается
			double delta_time = glfwGetTime() - tree_bolls[i].birth;
			glm::vec3 boll_centre = tree_bolls[i].centre + tree_bolls[i].direction * (float)delta_time * tree_bolls[i].speed;
			glm::vec3 player_position = getPosition();
			double distance = std::sqrt(std::pow(boll_centre.x - player_position.x, 2) +
				std::pow(boll_centre.y - player_position.y, 2) +
				std::pow(boll_centre.z - player_position.z, 2));
			if (distance < tree_bolls[i].radius) {
				// Удаление файербола
				std::swap(tree_bolls[i], tree_bolls[tree_bolls.size() - 1]);
				tree_bolls.pop_back();
				// Эффект для игрока.
				effect_start = glfwGetTime();
				break;
			}
		}
		
		// Обновление пола
		if (glfwGetTime() - prev_floor_build_time > 10) {
			//Удаление старого пола 
		    for (int i = 0; i < floor.size(); ++i) {
			    std::swap(floor[i], floor[floor.size() - 1]);
			    floor.pop_back();
		    }

			// Добавление нового пола
			floor.resize(floor.size() + 1);
			floor_maker.add_floor(getPosition(), 10, &floor[floor.size() - 1]);
			fl_vertex_buffer_data = {};
			fl_uv_buffer_data = {};
			for (Floor fl : floor) {
		    	fl_vertex_buffer_data.insert(fl_vertex_buffer_data.end(), fl.floor.begin(), fl.floor.end());
		        fl_uv_buffer_data.insert(fl_uv_buffer_data.end(), fl.uves.begin(), fl.uves.end());
		    }

			prev_floor_build_time = glfwGetTime();
			drawFloor(TextureID_floor, Texture_floor, floorvertexbuffer, flooruvbuffer, fl_vertex_buffer_data, fl_uv_buffer_data);
		}

		// Обновление неба
		if (glfwGetTime() - prev_sky_build_time > 20) {
			//Удаление старого неба 
			for (int i = 0; i < sky.size(); ++i) {
				std::swap(sky[i], sky[sky.size() - 1]);
				sky.pop_back();
			}

			// Добавление нового неба
			sky.resize(sky.size() + 1);
			sky_maker.add_sky(getPosition(), &sky[sky.size() - 1]);
			sky_vertex_buffer_data = {};
			sky_uv_buffer_data = {};
			for (Sky sk : sky) {
				sky_vertex_buffer_data.insert(sky_vertex_buffer_data.end(), sk.sky.begin(), sk.sky.end());
				sky_uv_buffer_data.insert(sky_uv_buffer_data.end(), sk.uves.begin(), sk.uves.end());
			}

			prev_sky_build_time = glfwGetTime();
			drawSky(TextureID_sky, Texture_sky, skyvertexbuffer, skyuvbuffer, sky_vertex_buffer_data, sky_uv_buffer_data);
		}

		// Действие эффекта попадания в игрока (в течении 5 секунд : 1я половина секунды - бездействие)
		if (glfwGetTime() - effect_start < 0.5 && alive == 1) {
			alive = 0;
			foreground.resize(foreground.size() + 1);
			foreground_maker.add_foreground(getPosition(), getDirection(), &foreground[foreground.size() - 1]);
		}

		if (glfwGetTime() - effect_start >= 0.5 && glfwGetTime() - effect_start < 5) {
			for (int i = 0; i < foreground.size(); ++i) {
				std::swap(foreground[i], foreground[foreground.size() - 1]);
				foreground.pop_back();
			}

			glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_TRUE);
			alive = 1;
		}

		// Возвращение в нормальный режим
		if (glfwGetTime() - effect_start >= 5) {
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		}

		// Добавление новых деревьев
		double end = glfwGetTime();
		if (end - start > interval) {
			trees.resize(trees.size() + 1);
			tree_maker.add_tree(&trees[trees.size() - 1]);
			interval = dis(gen) / 10000.0 * 5;
			start = glfwGetTime();
		}

		std::vector<GLfloat> g_vertex_buffer_data = {};
		std::vector<GLfloat> g_color_buffer_data = {};
		std::vector<GLfloat> g_uv_buffer_data = {};

		// Добавление всех объектов в буферы
		for (Fireboll fireboll : firebolls) {
			// Нужно учесть, что фаерболы двигаются
			Fireboll fireboll_copy = fireboll;
			fireboll_maker.find_position(&fireboll_copy, current_scene_speed);
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), fireboll_copy.boll.begin(), fireboll_copy.boll.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), fireboll_copy.boll_colors.begin(), fireboll_copy.boll_colors.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), fireboll_copy.uves.begin(), fireboll_copy.uves.end());
		}

		for (SpruceFireboll tree_boll : tree_bolls) {
			// Нужно учесть, что фаерболы двигаются
			SpruceFireboll tree_boll_copy = tree_boll;
			tree_boll_maker.find_position(&tree_boll_copy, current_spruce_speed);
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), tree_boll_copy.boll.begin(), tree_boll_copy.boll.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), tree_boll_copy.boll_colors.begin(), tree_boll_copy.boll_colors.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), tree_boll_copy.uves.begin(), tree_boll_copy.uves.end());
		}

		for (Spruce tree : trees) {
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), tree.tree.begin(), tree.tree.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), tree.tree_colors.begin(), tree.tree_colors.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), tree.uves.begin(), tree.uves.end());
		}

		for (Stump stump : stumps) {
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), stump.tree.begin(), stump.tree.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), stump.tree_colors.begin(), stump.tree_colors.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), stump.uves.begin(), stump.uves.end());
		}

		for (Foreground f : foreground) {
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), f.foreground.begin(), f.foreground.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), f.foreground_color.begin(), f.foreground_color.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), f.uves.begin(), f.uves.end());
		}

		for (ExplodedBall explosion : exploded_balls) {
			// Нужно учесть, что части фаербола двигаются
			ExplodedBall explosion_copy = explosion;
			fireboll_maker.find_explosion_coords(&explosion_copy, current_scene_speed);
			g_vertex_buffer_data.insert(g_vertex_buffer_data.end(), explosion_copy.boll.begin(), explosion_copy.boll.end());
			g_color_buffer_data.insert(g_color_buffer_data.end(), explosion_copy.boll_colors.begin(), explosion_copy.boll_colors.end());
			g_uv_buffer_data.insert(g_uv_buffer_data.end(), explosion_copy.uves.begin(), explosion_copy.uves.end());
		}

		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glBufferData(GL_ARRAY_BUFFER, g_vertex_buffer_data.size() * 4, g_vertex_buffer_data.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glBufferData(GL_ARRAY_BUFFER, g_color_buffer_data.size() * 4, g_color_buffer_data.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
		glBufferData(GL_ARRAY_BUFFER, g_uv_buffer_data.size() * 4, g_uv_buffer_data.data(), GL_STATIC_DRAW);

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID1);

		// Compute the MVP matrix from keyboard and mouse input
		computeMatricesFromInputs();

		glm::mat4 ProjectionMatrix = getProjectionMatrix();
		glm::mat4 ViewMatrix = getViewMatrix();
		glm::mat4 ModelMatrix = glm::mat4(1.0);
		glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, textures[current_texture_num]);
		glUniform1i(texture_ids[current_texture_num], 0);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			2,                                // size : U+V => 2
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// 3rd attribute : colors
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glVertexAttribPointer(
			2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw the triangle !
		glDrawArrays(GL_TRIANGLES, 0, g_vertex_buffer_data.size() / 3);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);

		// draw floor
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, floor_textures[current_floor_texture_num]);
		glUniform1i(TextureID_floor, 0);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, floorvertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, flooruvbuffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			2,                                // size : U+V => 2
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw the triangles 
		glDrawArrays(GL_TRIANGLES, 0, fl_vertex_buffer_data.size() / 3);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		// draw sky 
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, sky_textures[current_sky_texture_num]);
		glUniform1i(TextureID_sky, 0);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, skyvertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, skyuvbuffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			2,                                // size : U+V => 2
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw the triangles 
		glDrawArrays(GL_TRIANGLES, 0, sky_vertex_buffer_data.size() / 3);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		char text[256];
		sprintf(text, "SCORE:%d ", player_score);
		printText2D(text, 480, 550, 36);

		char text2[256];
		sprintf(text2, "YOU LOSE");
		if (alive == 0) {
			printText2D(text2, 100, 300, 70);
		}

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
		cleanupText2D;

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &colorbuffer);
	glDeleteBuffers(1, &uvbuffer);
	glDeleteProgram(programID1);
	glDeleteTextures(1, &Texture);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
