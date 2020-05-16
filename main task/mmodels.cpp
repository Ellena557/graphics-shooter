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
		tree_copy[i + 1] = tree_copy[i + 1] / scale + bias_y;
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

	res->turn_xy = 2 * PI * dis(gen) / 10000.0;
	res->turn_xz = 2 * PI * dis(gen) / 10000.0;
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

MakeFireboll::MakeFireboll()
{
	int N_PSI_part = N_PSI * 2 / 3;
	boll.resize(3 * 6 * N_PSI * N_PHI);
	boll_colors.resize(3 * 6 * N_PSI * N_PHI);
	uves.resize(2 * N_PSI * N_PHI * 6);

	glm::vec3 vertices[N_PSI][N_PHI];

	for (int j = 0; j < N_PSI; ++j) {
		for (int i = 0; i < N_PHI; ++i) {
			float phi = i * 2 * PI / (N_PHI - 1);
			float psi = j * PI / (N_PSI - 1);
			vertices[j][i] = glm::vec3(sin(psi) * cos(phi), sin(psi) * sin(phi), cos(psi));
		}
	}

	//add u, v, poles are with z coordinate
	glm::vec2 uv_coords[N_PSI][N_PHI];

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
}

void MakeFireboll::find_position(Fireboll* res)
{
	double delta_time = glfwGetTime() - res->birth;
	for (int i = 0; i < res->boll.size(); i += 3) {
		res->boll[i] += res->direction.x * delta_time * res->speed;
		res->boll[i + 1] += res->direction.y * delta_time * res->speed;
		res->boll[i + 2] += res->direction.z * delta_time * res->speed;
	}
	res->centre += res->direction * (float)delta_time * res->speed;
}

MakeSpruceFireboll::MakeSpruceFireboll()
{
	int N_PSI_part = N_PSI * 2 / 3;
	boll.resize(2 * N_PSI_part * N_PHI * 9);
	boll_colors.resize(2 * N_PSI_part * N_PHI * 9);

	glm::vec3 vertices[N_PSI][N_PHI];

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
		//boll_colors[i + 2] = 0.;
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

void MakeSpruceFireboll::find_position(SpruceFireboll* res)
{
	double delta_time = glfwGetTime() - res->birth;
	for (int i = 0; i < res->boll.size(); i += 3) {
		res->boll[i] += res->direction.x * delta_time * res->speed;
		res->boll[i + 1] += res->direction.y * delta_time * res->speed;
		res->boll[i + 2] += res->direction.z * delta_time * res->speed;
	}
	res->centre += res->direction * (float)delta_time * res->speed;
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