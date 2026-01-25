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

	ThreeVector operator*=(const double &alpha) {
		e[0] *= alpha;
		e[1] *= alpha;
		e[2] *= alpha;
		return *this;
	}

	ThreeVector operator+=(const ThreeVector &other) {
		e[0] += other.x();
		e[1] += other.y();
		e[2] += other.z();
		return *this;
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
	double attenuation; //Energy loss on reflect

	MaterialProperties(ThreeVector col, double diffuse, double attenuate) : colour(col), diffusivity(diffuse), attenuation(attenuate) {}
	MaterialProperties() : colour(ThreeVector(0, 0, 0)), diffusivity(0), attenuation(0.97) {}
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

		if (discriminant > 0) {
			double distance1 = (-half_b - std::sqrt(discriminant)) / a;
				if (distance1 > 0) {
					return distance1;
				}
		}
		return std::numeric_limits<double>::max();
	}
	ThreeVector normal(ThreeVector &position) const {
		return (position - centre) / radius;
	}
};

class Light {
public:
	ThreeVector position;
	double brightness;
	
	Light(ThreeVector pos, double bright) : position(pos), brightness(bright) {}
	Light() : position(ThreeVector()), brightness(0.8) {}

	ThreeVector skycolour(const ThreeVector &ray_dir) const {
		double a = 2 * (ray_dir.normalise().y() + 1);
		return ThreeVector(255, 255, 255) * (1 - a) + ThreeVector(180, 220, 255) * a;
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
	Light light;
	std::list<Sphere> spheres;

	Scene(Camera c, Light l, std::list<Sphere> s) : camera(c), light(l), spheres(s) {}
};

struct ray_payload {
	Sphere* sphere_pointer;
	double distance;
	ray_payload(Sphere* point, double dis) : sphere_pointer(point), distance(dis) {}
};

ray_payload get_payload(Ray& current_ray, const std::list<Sphere>& spheres){
	Sphere* current_object_pointer = nullptr;
	double min_distance = std::numeric_limits<double>::max();

	for (const Sphere& sphere : spheres) {
		double current_distance = sphere.intersect(current_ray);
		if (current_distance < min_distance) {
			min_distance = current_distance;
			const Sphere* current_object_pointer = &sphere;
		}
	}
	ray_payload payload(current_object_pointer, min_distance);
	return payload;
}

ThreeVector emit_ray(Ray& current_ray, const std::list<Sphere>& spheres, const Light& all_light, const int& remaining_bounces, const uint32_t& seed) {
	ThreeVector current_colour(0,0,0);
	double multiplier = 1;

	for (int i = 0; i < remaining_bounces; i++) {
		ray_payload payload = get_payload(current_ray, spheres);

		if (payload.sphere_pointer != nullptr) {
			ThreeVector intersect_pos = current_ray.origin + current_ray.direction*payload.distance;
			ThreeVector hitpoint_normal = payload.sphere_pointer->normal(intersect_pos);

			//General directional lighting
			ThreeVector light_to_sphere = (all_light.position - intersect_pos).normalise(); 
			double lightsource_reflection = all_light.brightness * hitpoint_normal.dot(light_to_sphere);
			
			//Attenuation of light for a more physical result// 
			/*
			ThreeVector light_to_sphere = (all_light.position - intersect_pos); 
			double lightsource_reflection = all_light.brightness * hitpoint_normal.dot(light_to_sphere)/light_to_sphere.length_squared();
			*/

			current_colour += payload.sphere_pointer->material.colour * (std::max(lightsource_reflection,0.2)) * multiplier;
			multiplier *= payload.sphere_pointer->material.attenuation;

			current_ray.origin = intersect_pos + hitpoint_normal*0.0001;
			ThreeVector new_normal = hitpoint_normal + ThreeVector(-0.5, 0.5, seed) * payload.sphere_pointer->material.diffusivity;
			current_ray.direction = current_ray.direction.reflect(new_normal);
		}
		else {
			current_colour += all_light.skycolour(current_ray.direction) * multiplier;
			break;
		}
	}
	return current_colour;
}

//convert linear index to 
std::tuple<int, int> to_array(const int i, const int width){
	int x = i % width;
	int y = (i - x) / width;
	return {x, y};
}

void render(std::vector<unsigned char>& pixels, const Scene& scene, const int& width, const int& height) {
	//Very expensive
	int bounces = 2;
	const int ray_count = 1;

	for (int i = 0; i < width * height; i++){
		const auto [x, y] = to_array(i, width);

		double u = (x - width / 2.0) / static_cast<double>(width);
		double v = -(y - height / 2.0) / static_cast<double>(height);

		ThreeVector pixel_colour(0, 0, 0);

		uint32_t seed = (y + 1) * (x + 1);
		double random_shift_x = 0.0;
		double random_shift_y = 0.0;

		for (int j = 0; j < ray_count; j++) {
			ThreeVector ray_direction = scene.camera.direction + scene.camera.right * (u + random_shift_x) + scene.camera.up * (v + random_shift_y);
			Ray ray1(scene.camera.position, ray_direction);
			pixel_colour += emit_ray(ray1, scene.spheres, scene.light, bounces, seed).clamp(0, 255);

			random_shift_x = (RandomDouble(seed) - 0.5) / width;
			random_shift_y = (RandomDouble(seed) - 0.5) / height;
		}
		pixel_colour = pixel_colour / ray_count;

		pixels[i * 3] = static_cast<int>(pixel_colour.x());
		pixels[i * 3 + 1] = static_cast<int>(pixel_colour.y());
		pixels[i * 3 + 2] = static_cast<int>(pixel_colour.z());
	}
}

int main() {
	const int width = 100;
	const int height = 100;
	const int array_size = width*height*3;

	std::string filename = "test";
	std::vector<unsigned char> pixels(array_size);

	Sphere object(ThreeVector(0, 0, 5), 3, MaterialProperties(ThreeVector(100, 100, 100), 0.0, 0.97));
	std::list<Sphere> spheres = { object };

	Light lightsource(ThreeVector(0, 0, 0), 200);
	Camera camera(ThreeVector(0, 0, 0), ThreeVector(0, 1, 0), ThreeVector(0, 0, 1));

	Scene scene(camera, lightsource, spheres);

	auto start = std::chrono::high_resolution_clock::now();
	render(pixels, scene, width, height);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::clog << "Render took " << duration.count() << " milliseconds." << std::endl;
	
	save_png(filename, width, height, pixels);
	return 0;
}
