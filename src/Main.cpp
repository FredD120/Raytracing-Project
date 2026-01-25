#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>
#include <cmath>
#include <list>
#include <limits>
#include <random>
#include <vector>
#include <chrono>
#include <tuple>
#include <filesystem>

//take in pixel colours as a vector of rgb values for a flattened array of size width * height
void save_png(const std::string filename, const int width, const int height, std::vector<unsigned char> pixels) {
	std::string name = filename + ".png";
	std::filesystem::path filepath = std::filesystem::current_path().parent_path() / "images" / name;

    stbi_write_png(filepath.string().c_str(), width, height, 3, pixels.data(), width * 3);
}

//returns a random 32-bit integer, using a fast RNG
static uint32_t PCGHash(uint32_t input) {
	uint32_t state = input * 747796405u + 2891336453u;
	uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
//modifies seed using RNG and returns a random float-32 between [0, 1]
static double RandomDouble(uint32_t& seed) {
	seed = PCGHash(seed);
	return ((double)seed / (double)std::numeric_limits<uint32_t>::max());
}

class ThreeVector {
public:
	double e[3];
	ThreeVector() : e{ 0,0,0 } {}

	ThreeVector(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

	ThreeVector(double& min, double& max, uint32_t& seed) { //fast rng
		e[0] = RandomDouble(seed) * (max - min) + min;
		e[1] = RandomDouble(seed) * (max - min) + min;
		e[2] = RandomDouble(seed) * (max - min) + min;
	}

	double x() const { return e[0]; }
	double y() const { return e[1]; }
	double z() const { return e[2]; }
	
	double dot(const ThreeVector &other) const {
		return e[0] * other.e[0] + e[1] * other.e[1] + e[2] * other.e[2];
	}

	ThreeVector cross(const ThreeVector &other) const {
		return ThreeVector(e[1] * other.e[2] - e[2] * other.e[1], e[2] * other.e[0] - e[0] * other.e[2], e[0] * other.e[1] - e[1] * other.e[0]);
	}

	void print() const {
		std::clog << e[0] << ' ' << e[1] << ' ' << e[2] << '\n';
	}

	ThreeVector normalise() const {
		double mag = magnitude();
		if (mag == 0){
			std::clog << "Magnitude zero vector cannot be normalised" << std::endl;
			return ThreeVector(e[0], e[1], e[2]);
		} else {
			return ThreeVector(e[0] / mag, e[1] / mag, e[2] / mag);
		}
	}

	double length_squared() const {
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	double magnitude() const {
		return std::sqrt(length_squared());
	}

	ThreeVector operator*=(const double& alpha) {
		e[0] *= alpha;
		e[1] *= alpha;
		e[2] *= alpha;
		return *this;
	}

	ThreeVector operator+=(const ThreeVector& other) {
		e[0] += other.x();
		e[1] += other.y();
		e[2] += other.z();
		return *this;
	}

	ThreeVector operator*=(const ThreeVector& other) {
		e[0] *= other.x();
		e[1] *= other.y();
		e[2] *= other.z();
		return *this;
	}

	ThreeVector operator*(const ThreeVector& other) const {
		return ThreeVector(e[0] * other.x(), e[1] * other.y(), e[2] * other.z());
	}

	ThreeVector operator*(const double alpha) const {
		return ThreeVector(e[0] * alpha, e[1] * alpha, e[2] * alpha);
	}

	ThreeVector operator/(const double alpha) const {
		return ThreeVector(e[0] / alpha, e[1] / alpha, e[2] / alpha);
	}

	ThreeVector operator-(const ThreeVector &other) const {
		return ThreeVector(e[0] - other.e[0], e[1] - other.e[1], e[2] - other.e[2]);
	}

	ThreeVector operator+(const ThreeVector &other) const {
		return ThreeVector(e[0] + other.e[0], e[1] + other.e[1], e[2] + other.e[2]);
	}

	ThreeVector invert() const {
		return ThreeVector(-e[0], -e[1], -e[2]);
	}

	ThreeVector clamp(double min, double max) {
		for (double &coord : e) {
			if (coord > max) {
				coord = max;
			}
			else if (coord < min) {
				coord = min;
			}
		}
		return *this;
	}

	ThreeVector reflect(ThreeVector& normal) {
		return *this - normal * ((2 * this->dot(normal) / normal.length_squared()));
	}
};

const ThreeVector white(1.0, 1.0, 1.0);
const ThreeVector blue(0.706, 0.863, 1.0);

ThreeVector skycolour(const ThreeVector &ray_dir) {
	double a = 2 * (ray_dir.normalise().y() + 1);
	return white * (1 - a) + blue * a;
}

class Ray {
public:
	ThreeVector origin;
	ThreeVector direction;

	Ray(ThreeVector o, ThreeVector d) : origin(o), direction(d) {}

	Ray() : origin(ThreeVector()), direction(ThreeVector()) {}
};

class MaterialProperties {
public:
	ThreeVector colour;
	double diffusivity; //Extent to which light is reflected deterministically
	double brightness; //Energy loss on reflect

	MaterialProperties(ThreeVector col, double diffuse, double bright) : colour(col), diffusivity(diffuse), brightness(bright) {}
	MaterialProperties() : colour(ThreeVector(0, 0, 0)), diffusivity(0), brightness(0.97) {}
};

class Sphere {
public:
	ThreeVector centre;
	double radius;
	MaterialProperties material;

	Sphere(ThreeVector cen, double rad, MaterialProperties mat) : centre(cen), radius(rad), material(mat) {}
	Sphere() : centre(ThreeVector()), radius(1), material(MaterialProperties()) {}

	double intersect(Ray current_ray) const {
		ThreeVector origin_to_center = current_ray.origin - centre;
		double a = current_ray.direction.length_squared();
		double half_b = origin_to_center.dot(current_ray.direction);
		double c = origin_to_center.dot(origin_to_center) - radius * radius;
		double discriminant = half_b * half_b - a * c;
		double sqrt_discrim = std::sqrt(discriminant);

		if (discriminant > 0) {
			double distance1 = (-half_b - sqrt_discrim) / a;
			if (distance1 > 0) {
				return distance1;
			}
			else{
				double distance2 = (-half_b + sqrt_discrim) / a;
				if (distance2 > 0) {
					return distance2;
				}
			}
		}
		return std::numeric_limits<double>::max();
	}
	ThreeVector normal(ThreeVector &position) const {
		return (position - centre) / radius;
	}
};

class Camera {
public:
	ThreeVector	position;
	ThreeVector up;
	ThreeVector direction;
	ThreeVector right;

	Camera(ThreeVector pos, ThreeVector u, ThreeVector dir) :position(pos), up(u), direction(dir), right(dir.cross(up)) {}
	Camera() :position(ThreeVector(0,0,0)), up(ThreeVector(0, 1, 0)), direction(ThreeVector(0, 0, 1)), right(ThreeVector(1, 0, 0)) {}
};

struct Scene { 
public:
	Camera camera;
	std::list<Sphere> spheres;

	Scene(Camera c, std::list<Sphere> s) : camera(c), spheres(s) {}
};

struct ray_payload {
	const Sphere* sphere_pointer;
	double distance;
	bool hit;
	ray_payload(const Sphere* point, double dis, bool h) : sphere_pointer(point), distance(dis), hit(h) {}
};

ray_payload get_payload(Ray& current_ray, const std::list<Sphere>& spheres){
	const Sphere* current_object_pointer = nullptr;
	double min_distance = std::numeric_limits<double>::max();

	for (const Sphere& sphere : spheres) {
		double current_distance = sphere.intersect(current_ray);
		if (current_distance < min_distance) {
			min_distance = current_distance;
			current_object_pointer = &sphere;
		}
	}

    bool hit = min_distance < std::numeric_limits<double>::max();
	ray_payload payload(current_object_pointer, min_distance, hit);
	return payload;
}

const float safety_factor = 0.0001;
ThreeVector emit_ray(Ray& current_ray, const std::list<Sphere>& spheres, const int& remaining_bounces, uint32_t& seed) {
	ThreeVector absorption(1.0, 1.0, 1.0);
	ThreeVector colour_buffer(0, 0, 0);

	for (int i = 0; i < remaining_bounces; i++) {
		ray_payload payload = get_payload(current_ray, spheres);

		if (payload.hit) {
			MaterialProperties material = payload.sphere_pointer->material;
			ThreeVector intersect_pos = current_ray.origin + current_ray.direction * payload.distance;
			ThreeVector hitpoint_normal = payload.sphere_pointer->normal(intersect_pos);

			current_ray.origin = intersect_pos + hitpoint_normal * safety_factor;
			ThreeVector new_normal = hitpoint_normal + ThreeVector(-0.5, 0.5, RandomDouble(seed)) * material.diffusivity;
			current_ray.direction = current_ray.direction.reflect(new_normal);

			absorption *= material.colour;
			colour_buffer += absorption * material.brightness;
		}
		else {
			colour_buffer += skycolour(current_ray.direction) * absorption;
			break;
		}
	}
	return colour_buffer;
}

//convert linear index to 
std::tuple<int, int> to_array(const int i, const int width){
	int x = i % width;
	int y = (i - x) / width;
	return {x, y};
}

void render(std::vector<unsigned char>& pixels, const Scene& scene, const int& width, const int& height) {
	//Very expensive
	int bounces = 1;
	const int ray_count = 1;

	for (int i = 0; i < width * height; i++) {
		const auto [x, y] = to_array(i, width);

		double u = (x - width / 2.0) / static_cast<double>(width);
		double v = -(y - height / 2.0) / static_cast<double>(height);

		ThreeVector pixel_colour(0, 0, 0);

		uint32_t seed = static_cast<uint32_t>(i);
		double random_shift_x = 0.0;
		double random_shift_y = 0.0;

		for (int j = 0; j < ray_count; j++) {
			ThreeVector ray_direction = scene.camera.direction + scene.camera.right * (u + random_shift_x) + scene.camera.up * (v + random_shift_y);
			Ray ray1(scene.camera.position, ray_direction);
			pixel_colour += emit_ray(ray1, scene.spheres, bounces, seed);

			random_shift_x = (RandomDouble(seed) - 0.5) / width;
			random_shift_y = (RandomDouble(seed) - 0.5) / height;
		}
		pixel_colour = pixel_colour / ray_count;

		pixels[i * 3] = static_cast<int>(pixel_colour.x() * 255);
		pixels[i * 3 + 1] = static_cast<int>(pixel_colour.y() * 255);
		pixels[i * 3 + 2] = static_cast<int>(pixel_colour.z() * 255);
	}
}

int main() {
	const int width = 500;
	const int height = 500;
	const int array_size = width*height*3;

	std::string filename = "test";
	std::vector<unsigned char> pixels(array_size);

	Sphere object(ThreeVector(0, 0, 5), 1, MaterialProperties(ThreeVector(0.2, 0.8, 0.5), 1.0, 0.97));
	std::list<Sphere> spheres = { object };
	Camera camera(ThreeVector(0, 0, 0), ThreeVector(0, 1, 0), ThreeVector(0, 0, 1));

	Scene scene(camera, spheres);

	auto start = std::chrono::high_resolution_clock::now();
	render(pixels, scene, width, height);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::clog << "Render took " << duration.count() << " milliseconds." << std::endl;
	
	save_png(filename, width, height, pixels);
	return 0;
}
