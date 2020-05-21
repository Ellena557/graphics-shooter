#include "models.hpp"
#include <math.h>       

void make_tree(Tree* res, const std::vector<GLfloat>& tree, const std::vector<GLfloat>& tree_colors,
	double turn_xy, double turn_xz, double scale, double bias_x, double bias_y, double bias_z)
{
	// Новое дерево
	std::vector<GLfloat> tree_copy = tree;

	// Поворот в плоскости X-Y
	for (int i = 0; i < tree_copy.size(); i += 3) {
		float a = tree_copy[i];
		float b = tree_copy[i + 1];
		tree_copy[i] = a * cos(turn_xy) - b * sin(turn_xy);
		tree_copy[i + 1] = b * cos(turn_xy) + a * sin(turn_xy);
	}

	// Поворот в плоскости X-Z
	for (int i = 0; i < tree_copy.size(); i += 3) {
		float a = tree_copy[i];
		float b = tree_copy[i + 2];
		tree_copy[i] = a * cos(turn_xz) - b * sin(turn_xz);
		tree_copy[i + 2] = b * cos(turn_xz) + a * sin(turn_xz);
	}

	// Смещение и масштабирование.
	for (int i = 0; i < tree_copy.size(); i += 3) {
		tree_copy[i] = tree_copy[i] / scale + bias_x;
		tree_copy[i + 1] = tree_copy[i + 1] / scale; // +bias_y;
		tree_copy[i + 2] = tree_copy[i + 2] / scale + bias_z;
	}

	res->tree = tree_copy;
	res->tree_colors = tree_colors;
}


void MakeSpruce::add_tree(Spruce* res)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 10000);

	res->turn_xy = 2 * PI * dis(gen) / 150000.0;
	res->turn_xz = 2 * PI * dis(gen) / 150000.0;
	res->scale = dis(gen) / 2000.0 + 5;
	res->bias_x = dis(gen) / 500.0 - 10;
	res->bias_y = dis(gen) / 500.0 - 10;
	res->bias_z = dis(gen) / 500.0 - 10;

	make_tree(res, tree, tree_colors, res->turn_xy, res->turn_xz, res->scale, res->bias_x, res->bias_y, res->bias_z);

	res->radius = 0.3 / res->scale;
	res->centre = { res->tree[108], res->tree[109], res->tree[110] };
}

void MakeStump::add_tree(Stump* res, const Spruce& spruce)
{
	// Все параметры как у соответствующего дерева.
	make_tree(res, stump, stump_colors, spruce.turn_xy, spruce.turn_xz, spruce.scale, spruce.bias_x, spruce.bias_y, spruce.bias_z);
	res->birth = glfwGetTime();
}

//MakeFireboll::MakeFireboll(int N_PHI, int N_PSI)
MakeFireboll::MakeFireboll(int detals)
{
	int N_PHI = detals;
	int N_PSI = detals;
	int N_PSI_part = N_PSI * 2 / 3;
	boll.resize(3 * 6 * N_PSI * N_PHI);
	boll_colors.resize(3 * 6 * N_PSI * N_PHI);
	uves.resize(2 * N_PSI * N_PHI * 6);
	normals.resize(3 * N_PSI * N_PHI * 6);

	glm::vec3 vertices[100][100];

	for (int j = 0; j < N_PSI; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			float phi = i * 2 * PI / (N_PHI - 1);
			float psi = j * PI / (N_PSI - 1);
			vertices[j][i] = glm::vec3(sin(psi) * cos(phi), sin(psi) * sin(phi), cos(psi));
		}
	}

	//add u, v, poles are with z coordinate
	glm::vec2 uv_coords[100][100];

	for (int j = 1; j < N_PSI_part; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			float phi = i * 2 * PI / N_PHI;
			float psi = j * PI / N_PSI;

			float u = phi / (2 * PI) + 0.5;
			float v = psi / PI + 0.5;
			uv_coords[j][i] = glm::vec2(u, v);
		}
	}

	int k = 0;
	int n = 0;

	// add all sphere triangles
	for (int j = 0; j < N_PSI - 1; ++j) {
		for (int i = 0; i < N_PHI; ++i) {

			boll[k++] = vertices[j][i].x;
			boll[k++] = vertices[j][i].y;
			boll[k++] = vertices[j][i].z;

			boll[k++] = vertices[j + 1][i].x;
			boll[k++] = vertices[j + 1][i].y;
			boll[k++] = vertices[j + 1][i].z;

			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].z;

			boll[k++] = vertices[j][i].x;
			boll[k++] = vertices[j][i].y;
			boll[k++] = vertices[j][i].z;

			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].z;

			boll[k++] = vertices[j][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j][(i + 1) % N_PHI].z;
		}
	}

	for (int i = 0; i < k; ++i) {
		boll_colors[i] = 0;
	}

	for (int i = 0; i < k; ++i) {
		normals[i] = 0;
	}

	// add texture  UV coordinates
	for (int j = 0; j < N_PSI - 1; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			uves[n++] = (atan2(vertices[j][i].y, vertices[j][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][i].y, vertices[j + 1][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][(i + 1) % N_PHI].y, vertices[j + 1][(i + 1) % N_PHI].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][(i + 1) % N_PHI].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j][i].y, vertices[j][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][(i + 1) % N_PHI].y, vertices[j + 1][(i + 1) % N_PHI].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][(i + 1) % N_PHI].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j][(i + 1) % N_PHI].y, vertices[j][(i + 1) % N_PHI].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][(i + 1) % N_PHI].z) / PI + 0.5;

			if (uves[n - 2] < uves[n - 6]) {
				uves[n - 8] = 1;
				uves[n - 4] = 1;
				uves[n - 2] = 1;
			}
		}
	}

}

MakeFireboll::MakeFireboll(std::vector<GLfloat> boll, std::vector<GLfloat> boll_colors, std::vector<GLfloat> uves, std::vector<GLfloat> normals) {
	
}

void MakeFireboll::add_boll(glm::vec3 position, glm::vec3 direction, int scale, Fireboll* res)
{
	// Нормализация направления.
	double boll_direction_length = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
	direction.x /= boll_direction_length;
	direction.y /= boll_direction_length;
	direction.z /= boll_direction_length;

	std::vector<GLfloat> boll_copy = boll;
	float a = direction.x;
	float b = direction.y;
	float c = direction.z;
	float d = sqrt(a * a + c * c);

	for (int i = 0; i < boll_copy.size(); ++i) {
		boll_copy[i] /= scale;
	}

	glm::vec3 z = direction;
	glm::vec3 x = { c / d, 0., -a / d };
	d = sqrt(b * b * x.z * x.z + (c * x.x - a * x.z) * (c * x.x - a * x.z) + b * b * x.x * x.x);
	glm::vec3 y = { b * x.z / d, (c * x.x - a * x.z) / d, -b * x.x / d };
	// Поворот
	for (int i = 0; i < boll_copy.size(); i += 3) {
		float p = boll_copy[i];
		float q = boll_copy[i + 1];
		float r = boll_copy[i + 2];
		boll_copy[i] = p * x.x + q * y.x + r * z.x;
		boll_copy[i + 1] = p * x.y + q * y.y + r * z.y;
		boll_copy[i + 2] = p * x.z + q * y.z + r * z.z;
	}

	// Смещение
	for (int i = 0; i < boll_copy.size(); i += 3) {
		boll_copy[i] += position.x;
		boll_copy[i + 1] += position.y;
		boll_copy[i + 2] += position.z;
	}

	res->boll = boll_copy;
	res->boll_colors = boll_colors;
	res->uves = uves;
	res->direction = direction;
	res->radius = 1.0 / scale;
	res->centre = position;
	res->birth = glfwGetTime();
	res->isExploded = 0;
	res->normals = normals;
}

void MakeFireboll::find_position(Fireboll* res, float speed)
{
	double delta_time = glfwGetTime() - res->birth;
	for (int i = 0; i < res->boll.size(); i += 3) {
		res->boll[i] += res->direction.x * delta_time * speed;
		res->boll[i + 1] += res->direction.y * delta_time * speed;
		res->boll[i + 2] += res->direction.z * delta_time * speed;
	}
	res->centre += res->direction * (float)delta_time * speed;
}

MakeSpruceFireboll::MakeSpruceFireboll(int N_PHI, int N_PSI)
{
	int N_PSI_part = N_PSI * 2 / 3;
	boll.resize(2 * N_PSI_part * N_PHI * 9);
	boll_colors.resize(2 * N_PSI_part * N_PHI * 9);

	glm::vec3 vertices[100][100];

	for (int j = 1; j < N_PSI_part; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			float phi = i * 2 * PI / N_PHI;
			float psi = j * PI / N_PSI;
			vertices[j][i] = glm::vec3(sin(psi) * cos(phi), sin(psi) * sin(phi), cos(psi));
		}
	}

	glm::vec3 north_pole = vec3(0., 0., 1.);
	glm::vec3 south_pole = vec3(0., 0., -2);

	int k = 0;

	//add triangles for south pole
	for (int i = 0; i < N_PHI; ++i) {
		boll[k++] = north_pole.x;
		boll[k++] = north_pole.y;
		boll[k++] = north_pole.z;
		boll_colors[k - 2] = north_pole.z;

		boll[k++] = vertices[1][i].x;
		boll[k++] = vertices[1][i].y;
		boll[k++] = vertices[1][i].z;
		boll_colors[k - 2] = vertices[1][i].z;

		boll[k++] = vertices[1][(i + 1) % N_PHI].x;
		boll[k++] = vertices[1][(i + 1) % N_PHI].y;
		boll[k++] = vertices[1][(i + 1) % N_PHI].z;
		boll_colors[k - 2] = vertices[1][(i + 1) % N_PHI].z;
	}

	//add triangles for north pole
	for (int i = 0; i < N_PHI; ++i) {
		boll[k++] = south_pole.x;
		boll[k++] = south_pole.y;
		boll[k++] = south_pole.z;
		boll_colors[k - 3] = 0.2;

		boll[k++] = vertices[N_PSI_part - 1][i].x;
		boll[k++] = vertices[N_PSI_part - 1][i].y;
		boll[k++] = vertices[N_PSI_part - 1][i].z;
		boll_colors[k - 2] = vertices[N_PSI_part - 1][i].z;


		boll[k++] = vertices[N_PSI_part - 1][(i + 1) % N_PHI].x;
		boll[k++] = vertices[N_PSI_part - 1][(i + 1) % N_PHI].y;
		boll[k++] = vertices[N_PSI_part - 1][(i + 1) % N_PHI].z;
		boll_colors[k - 2] = vertices[N_PSI_part - 1][(i + 1) % N_PHI].z;
	}

	//add all sphere triangles
	for (int j = 1; j < N_PSI_part - 1; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			boll[k++] = vertices[j][i].x;
			boll[k++] = vertices[j][i].y;
			boll[k++] = vertices[j][i].z;
			boll_colors[k - 2] = vertices[j][i].z;

			boll[k++] = vertices[j + 1][i].x;
			boll[k++] = vertices[j + 1][i].y;
			boll[k++] = vertices[j + 1][i].z;
			boll_colors[k - 2] = vertices[j + 1][i].z;

			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].z;
			boll_colors[k - 2] = vertices[j + 1][(i + 1) % N_PHI].z;

			//////

			boll[k++] = vertices[j][i].x;
			boll[k++] = vertices[j][i].y;
			boll[k++] = vertices[j][i].z;
			boll_colors[k - 2] = vertices[j][i].z;

			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j + 1][(i + 1) % N_PHI].z;
			boll_colors[k - 2] = vertices[j + 1][(i + 1) % N_PHI].z;

			boll[k++] = vertices[j][(i + 1) % N_PHI].x;
			boll[k++] = vertices[j][(i + 1) % N_PHI].y;
			boll[k++] = vertices[j][(i + 1) % N_PHI].z;
			boll_colors[k - 2] = vertices[j][(i + 1) % N_PHI].z;
		}
	}

	for (int i = 0; i < k; ++i) {
		if (boll_colors[i] < 0)
			boll_colors[i] = -boll_colors[i];
	}

	// Установка красного и синего цветов. Зелёный уже установлен на равномерное убывание от полюсов.
	for (int i = 0; i < k; i += 3) {
		boll_colors[i] = 1. - boll_colors[i]; // обычно до этого был 0.
	}

	for (int i = 0; i < k; ++i) {
		boll[i] /= scale;
	}
}
void MakeSpruceFireboll::add_boll(glm::vec3 position, glm::vec3 direction, SpruceFireboll* res)
{
	// Нормализация направления.
	double boll_direction_length = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
	direction.x /= boll_direction_length;
	direction.y /= boll_direction_length;
	direction.z /= boll_direction_length;

	std::vector<GLfloat> boll_copy = boll;
	float a = direction.x;
	float b = direction.y;
	float c = direction.z;
	float d = sqrt(a * a + c * c);

	glm::vec3 z = direction;
	glm::vec3 x = { c / d, 0., -a / d };
	d = sqrt(b * b * x.z * x.z + (c * x.x - a * x.z) * (c * x.x - a * x.z) + b * b * x.x * x.x);
	glm::vec3 y = { b * x.z / d, (c * x.x - a * x.z) / d, -b * x.x / d };
	//Поворот
	for (int i = 0; i < boll_copy.size(); i += 3) {
		float p = boll_copy[i];
		float q = boll_copy[i + 1];
		float r = boll_copy[i + 2];
		boll_copy[i] = p * x.x + q * y.x + r * z.x;
		boll_copy[i + 1] = p * x.y + q * y.y + r * z.y;
		boll_copy[i + 2] = p * x.z + q * y.z + r * z.z;
	}

	//Смещение
	for (int i = 0; i < boll_copy.size(); i += 3) {
		boll_copy[i] += position.x;
		boll_copy[i + 1] += position.y;
		boll_copy[i + 2] += position.z;
	}

	for (int i = 0; i < boll_copy.size(); i += 3) {
		res->uves.emplace_back(0.0f);
		res->uves.emplace_back(0.0f);
	}

	res->boll = boll_copy;
	res->boll_colors = boll_colors;
	res->direction = direction;
	res->radius = 1.0 / scale;
	res->centre = position;
	res->birth = glfwGetTime();
}

void MakeSpruceFireboll::find_position(SpruceFireboll* res, float speed)
{
	double delta_time = glfwGetTime() - res->birth;
	for (int i = 0; i < res->boll.size(); i += 3) {
		res->boll[i] += res->direction.x * delta_time * speed;
		res->boll[i + 1] += res->direction.y * delta_time * speed;
		res->boll[i + 2] += res->direction.z * delta_time * speed;
	}
	res->centre += res->direction * (float)delta_time * speed;
}

void MakeForeground::add_foreground(glm::vec3 position, glm::vec3 direction, Foreground* res)
{
	float a = direction.x;
	float b = direction.y;
	float c = direction.z;

	float xx = position.x;
	float yy = position.y;
	float zz = position.z;

	float d = -a * xx - b * yy - c * zz;

	// вершины квадрата
	std::vector<GLfloat> foreground = {
		xx - 1000.0f + 1.1f * a, yy - 1000.0f + 1.1f * b, zz + (1000.0f * a + 1000.0f * b) / c + 1.1f * c,
		xx + 1000.0f + 1.1f * a, yy + 1000.0f + 1.1f * b, zz - (1000.0f * a + 1000.0f * b) / c + 1.1f * c,
		xx - 1000.0f + 1.1f * a, yy + 1000.0f + 1.1f * b, zz + (1000.0f * a - 1000.0f * b) / c + 1.1f * c,

		xx - 1000.0f + 1.1f * a, yy - 1000.0f + 1.1f * b, zz + (1000.0f * a + 1000.0f * b) / c + 1.1f * c,
		xx + 1000.0f + 1.1f * a, yy + 1000.0f + 1.1f * b, zz - (1000.0f * a + 1000.0f * b) / c + 1.1f * c,
		xx + 1000.0f + 1.1f * a, yy - 1000.0f + 1.1f * b, zz + (- 1000.0f * a + 1000.0f * b) / c + 1.1f * c,
	};

	std::vector<GLfloat> foreground_color = {
		0.2f, 0.59f, 0.21f,
		0.2f, 0.59f, 0.21f,
		0.2f, 0.59f, 0.21f,

		0.2f, 0.59f, 0.21f,
		0.2f, 0.59f, 0.21f,
		0.2f, 0.59f, 0.21f,
	};

	std::vector<GLfloat> uves = {
		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f,

		0.0f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f
	};

	res->foreground = foreground;
	res->foreground_color = foreground_color;
	res->uves = uves;
}


void MakeFloor::add_floor(glm::vec3 position, int floor_range, Floor* res)
{
	float xx = position.x;
	float yy = position.y;
	float zz = position.z;
	int scale = 1;
	const int n_squares = 80*80;
	int f_range = 80;
	std::vector<GLfloat> floor(2 * (f_range+1) * 2 * (f_range+1)*3*3*2);
	std::vector<GLfloat> colors(2 * (f_range + 1) * 2 * (f_range + 1) * 3*3*2);
	std::vector<GLfloat> uv_coords(2 * (f_range + 1) * 2 * (f_range + 1) * 3*2*2);
	int counter = 0;
	int uv_counter = 0;
	// вершины квадратов
	for (int i = -f_range; i <= f_range; ++i) {
		for (int j = -f_range; j <= f_range; ++j) {
			//x
			floor[counter] = xx + scale * (i);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i+j);
			colors[counter] = 1.0f;
			counter++;

			//x
			floor[counter] = xx + scale * (i + 1);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i + j);
			colors[counter] = 1.0f;
			counter++;

			//x
			floor[counter] = xx + scale * (i + 1);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i + 1 + j + 1);
			colors[counter] = 1.0f;
			counter++;

			//x
			floor[counter] = xx + scale * (i);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i + j);
			colors[counter] = 1.0f;
			counter++;

			//x
			floor[counter] = xx + scale * (i);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i + j + 1);
			colors[counter] = 1.0f;
			counter++;

			//x
			floor[counter] = xx + scale * (i + 1);
			colors[counter] = 0.0f;
			counter++;

			//y
			floor[counter] = -1.0;
			colors[counter] = 0.0f;
			counter++;

			//z
			floor[counter] = zz + scale * (i + 1 + j + 1);
			colors[counter] = 1.0f;
			counter++;


			uv_coords[uv_counter] = 0.0f;
			uv_counter++;
			uv_coords[uv_counter] = 0.0f;
			uv_counter++;

			uv_coords[uv_counter] = 0.5f;
			uv_counter++;
			uv_coords[uv_counter] = 0.0f;
			uv_counter++;

			uv_coords[uv_counter] = 1.0f;
			uv_counter++;
			uv_coords[uv_counter] = 1.0f;
			uv_counter++;

			uv_coords[uv_counter] = 0.0f;
			uv_counter++;
			uv_coords[uv_counter] = 0.0f;
			uv_counter++;

			uv_coords[uv_counter] = 0.0f;
			uv_counter++;
			uv_coords[uv_counter] = 0.5f;
			uv_counter++;

			uv_coords[uv_counter] = 1.0f;
			uv_counter++;
			uv_coords[uv_counter] = 1.0f;
			uv_counter++;
		}
	}

	res->floor = floor;
	res->colors = colors;
	res->uves = uv_coords;
}


void MakeSky::add_sky(glm::vec3 position, Sky* res)
{
	sky.resize(3 * 6 * N_PSI_big * N_PHI_big);
	uves.resize(2 * N_PSI_big * N_PHI_big * 6);

	glm::vec3 vertices[N_PSI_big][N_PHI_big];

	for (int j = 0; j < N_PSI_big; ++j) {
		for (int i = 0; i < N_PHI_big; ++i) {
			float phi = i * 2 * PI / (N_PHI_big - 1);
			float psi = j * PI / (N_PSI_big - 1);
			vertices[j][i] = glm::vec3(sin(psi) * cos(phi), sin(psi) * sin(phi), cos(psi));
		}
	}

	int k = 0;
	int n = 0;

	// add all sphere triangles
	for (int j = 0; j < N_PSI_big - 1; ++j) {
		for (int i = 0; i < N_PHI_big; ++i) {
			sky[k++] = vertices[j][i].x;
			sky[k++] = vertices[j][i].y;
			sky[k++] = vertices[j][i].z;

			sky[k++] = vertices[j + 1][i].x;
			sky[k++] = vertices[j + 1][i].y;
			sky[k++] = vertices[j + 1][i].z;

			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].x;
			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].y;
			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].z;

			sky[k++] = vertices[j][i].x;
			sky[k++] = vertices[j][i].y;
			sky[k++] = vertices[j][i].z;

			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].x;
			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].y;
			sky[k++] = vertices[j + 1][(i + 1) % N_PHI_big].z;

			sky[k++] = vertices[j][(i + 1) % N_PHI_big].x;
			sky[k++] = vertices[j][(i + 1) % N_PHI_big].y;
			sky[k++] = vertices[j][(i + 1) % N_PHI_big].z;
		}
	}

	// add texture  UV coordinates
	for (int j = 0; j < N_PSI_big - 1; ++j) {
		for (int i = 0; i < N_PHI_big; ++i) {
			uves[n++] = (atan2(vertices[j][i].y, vertices[j][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][i].y, vertices[j + 1][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][(i + 1) % N_PHI_big].y, vertices[j + 1][(i + 1) % N_PHI_big].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][(i + 1) % N_PHI_big].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j][i].y, vertices[j][i].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][i].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j + 1][(i + 1) % N_PHI_big].y, vertices[j + 1][(i + 1) % N_PHI_big].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j + 1][(i + 1) % N_PHI_big].z) / PI + 0.5;

			uves[n++] = (atan2(vertices[j][(i + 1) % N_PHI_big].y, vertices[j][(i + 1) % N_PHI_big].x) / (2 * PI) + 0.5);
			uves[n++] = -asin(vertices[j][(i + 1) % N_PHI_big].z) / PI + 0.5;

			if (uves[n - 2] < uves[n - 6]) {
				uves[n - 8] = 1;
				uves[n - 4] = 1;
				uves[n - 2] = 1;
			}
		}
	}

	for (int i = 0; i < sky.size(); ++i) {
		sky[i] *= 70;
	}

	//Поворот
	for (int i = 0; i < sky.size(); i += 3) {
		float p = sky[i];
		float q = sky[i + 1];
		sky[i] = q;
		sky[i + 1] = p;
	}

	res->sky = sky;
	res->uves = uves;
}


void MakeFireboll::add_exploded_ball(ExplodedBall* res, Fireboll* fireball, double time, float speed, float normals_speed)
{
	res->speed = fireball->speed;
	int n = fireball->boll.size();
	std::vector<GLfloat> normals = fireball->normals;

	// direction of the normal
	float xx;
	float yy;
	float zz;
	int counter = 0;

	for (int i = 0; i < n; i += 9) {
		float x1 = fireball->boll[i];
		float y1 = fireball->boll[i + 1];
		float z1 = fireball->boll[i + 2];
		float x2 = fireball->boll[i + 3];
		float y2 = fireball->boll[i + 4];
		float z2 = fireball->boll[i + 5];
		float x3 = fireball->boll[i + 6];
		float y3 = fireball->boll[i + 7];
		float z3 = fireball->boll[i + 8];

		xx = (y3 - y1) * (z2 - z1) - (y2 - y1) * (z3 - z1);
		yy = -((x3 - x1) * (z2 - z1) - (x2 - x1) * (z3 - z1));
		zz = (x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1);
		
		normals[counter] = xx;
		normals[counter + 1] = yy;
		normals[counter + 2] = zz;

		float ddd = -xx * x1 - yy * y1 - zz * z1;

		glm::vec3 pos_centre = fireball->centre;

		if (pos_centre.x * xx + pos_centre.y * yy + pos_centre.z * zz + ddd > 0) {
			normals[counter] = -xx;
			normals[counter + 1] = -yy;
			normals[counter + 2] = -zz;
		}

		counter += 3;
	}

	for (int i = 0; i < normals.size(); i += 1) {
		//normals[i] *= 1000;
		normals[i] *= 50 * normals_speed * normals_speed;
	}

	double delta_time = time - fireball->birth;
	res->boll = fireball->boll;

	for (int i = 0; i < res->boll.size(); i += 3) {
		res->boll[i] += fireball->direction.x * delta_time * speed;
		res->boll[i + 1] += fireball->direction.y * delta_time * speed;
		res->boll[i + 2] += fireball->direction.z * delta_time * speed;
	}

	res->boll_colors = fireball->boll_colors;

	for (int i = 0; i < res->boll_colors.size(); i += 9) {
		res->boll_colors[i] = 1.0f;
		res->boll_colors[i + 1] = 0.01f;
		res->boll_colors[i + 2] = 0.01f;
		res->boll_colors[i + 3] = 0.98f;
		res->boll_colors[i + 4] = 0.79f;
		res->boll_colors[i + 5] = 0.02f;
		res->boll_colors[i + 6] = 0.98f;
		res->boll_colors[i + 7] = 0.32f;
		res->boll_colors[i + 8] = 0.02f;
	}

	res->uves = fireball->uves;
	res->normals = normals;
	res->birth = fireball->explosion_time;
}


void MakeFireboll::find_explosion_coords(ExplodedBall* res, float speed)
{
	double delta_time = glfwGetTime() - res->birth;

	int n = res->boll.size();
	int counter = 0;
	float xx;
	float yy;
	float zz;
	for (int i = 0; i < n; i += 9) {
		xx = res->normals[counter];
		yy = res->normals[counter + 1];
		zz = res->normals[counter + 2];
		counter += 3;

		res->boll[i] += xx * delta_time * speed;
		res->boll[i + 1] += yy * delta_time * speed;
		res->boll[i + 2] += zz * delta_time * speed;
		res->boll[i + 3] += xx * delta_time * speed;
		res->boll[i + 4] += yy * delta_time * speed;
		res->boll[i + 5] += zz * delta_time * speed;
		res->boll[i + 6] += xx * delta_time * speed;
		res->boll[i + 7] += yy * delta_time * speed;
		res->boll[i + 8] += zz * delta_time * speed;
	}
}
