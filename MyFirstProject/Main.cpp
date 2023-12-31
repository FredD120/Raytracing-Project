#include <iostream>
#include <cmath>
#include <list>
#include <limits>
#include <random>
#include <vector>
#include <chrono>

static uint32_t PCGHash(uint32_t input) { //faster rng
	uint32_t state = input * 747796405u + 2891336453u;
	uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
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
	
	double dot(const ThreeVector &other) {
		return e[0] * other.e[0] + e[1] * other.e[1] + e[2] * other.e[2];
	}

	ThreeVector cross(const ThreeVector &other) {
		return ThreeVector(e[1] * other.e[2] - e[2] * other.e[1], e[2] * other.e[0] - e[0] * other.e[2], e[0] * other.e[1] - e[1] * other.e[0]);
	}

	void print() {
		std::clog << e[0] << ' ' << e[1] << ' ' << e[2] << '\n';
	}

	ThreeVector normalise() {
		double mag = magnitude();
		if (mag == 0){
			std::clog << "Magnitude zero vector cannot be normalised" << std::endl;
			return *this;
		} else {
			return ThreeVector(e[0] / mag, e[1] / mag, e[2] / mag);
		}
	}

	double length_squared() {
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	double magnitude() {
		return std::sqrt(length_squared());
	}

	ThreeVector& operator*=(double &alpha) {
		e[0] *= alpha;
		e[1] *= alpha;
		e[2] *= alpha;
		return *this;
	}

	ThreeVector& operator+=(const ThreeVector &other) {
		e[0] += other.x();
		e[1] += other.y();
		e[2] += other.z();
		return *this;
	}

	ThreeVector operator*(double alpha) {
		return ThreeVector(e[0] * alpha, e[1] * alpha, e[2] * alpha);
	}

	ThreeVector operator/(double alpha) {
		return ThreeVector(e[0] / alpha, e[1] / alpha, e[2] / alpha);
	}

	ThreeVector operator-(const ThreeVector &other) {
		return ThreeVector(e[0] - other.e[0], e[1] - other.e[1], e[2] - other.e[2]);
	}

	ThreeVector operator+(const ThreeVector &other) {
		return ThreeVector(e[0] + other.e[0], e[1] + other.e[1], e[2] + other.e[2]);
	}

	ThreeVector invert() {
		return ThreeVector(-e[0], -e[1], -e[2]);
	}

	ThreeVector clamp(double min,double max) {
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
		return *this - normal*((2 * this->dot(normal) / normal.length_squared()));
	}
};

class Ray {
public:
	ThreeVector origin;
	ThreeVector direction;

	Ray(ThreeVector& o, ThreeVector& d) : origin(o), direction(d) {}

	Ray() : origin(ThreeVector()), direction(ThreeVector()) {}
};

class MaterialProperties {
public:
	ThreeVector colour;
	double diffusivity; //Extent to which light is reflected deterministically
	double attenuation; //Energy loss on reflect

	MaterialProperties(ThreeVector col, double diffuse, double attenuate) : colour(col), diffusivity(diffuse), attenuation(attenuate) {}
	MaterialProperties() : colour(ThreeVector(0,0,0)), diffusivity(0), attenuation(0.97) {}
};

class Sphere {
public:
	ThreeVector centre;
	double radius;
	MaterialProperties material;

	Sphere(ThreeVector cen, double rad, MaterialProperties mat) : centre(cen), radius(rad), material(mat){}
	Sphere() : centre(ThreeVector()), radius(1), material(MaterialProperties()) {}

	double intersect(Ray current_ray) {
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
			/*
			double distance1 = (-b - std::sqrt(discriminant)) / (2.0 * a); //alternate way to calculate using both solutions
			double distance2 = (-b + std::sqrt(discriminant)) / (2.0 * a);
			if (distance1 > 0 && distance2 > distance1) {
				return distance2;
			}
			else if (distance1 > 0) {
				return distance1;
			}
			*/
		}
		return std::numeric_limits<double>::max();
	}
	ThreeVector normal(ThreeVector &position) {
		return (position - centre)/radius;
	}
};

class Light {
public:
	
	ThreeVector position;
	double brightness;
	
	Light(ThreeVector pos, double bright) : position(pos), brightness(bright) {}
	Light() : position(ThreeVector()), brightness(0.8) {}

	ThreeVector skycolour(ThreeVector &ray_dir) {
		double a = 2*(ray_dir.normalise().y() + 1);
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

struct ray_payload {
	Sphere* sphere_pointer;
	double distance;
	ray_payload(Sphere* point, double dis) : sphere_pointer(point), distance(dis) {}
};

ray_payload get_payload(Ray& current_ray,std::list<Sphere>& spheres){
	Sphere* current_object_pointer = nullptr;
	double min_distance = std::numeric_limits<double>::max();
	for (Sphere& sphere : spheres) {
		double current_distance = sphere.intersect(current_ray);
		if (current_distance < min_distance) {
			min_distance = current_distance;
			current_object_pointer = &sphere;
		}
	}
	ray_payload payload(current_object_pointer, min_distance);
	return payload;
}

ThreeVector emit_ray(Ray &current_ray, std::list<Sphere>& spheres, Light &all_light, int &remaining_bounces, uint32_t &seed) {
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

			
			current_colour += payload.sphere_pointer->material.colour*(std::max(lightsource_reflection,0.2))*multiplier;
			multiplier *= payload.sphere_pointer->material.attenuation;

			current_ray.origin = intersect_pos + hitpoint_normal*0.0001;
			ThreeVector new_normal = hitpoint_normal + ThreeVector(-0.5,0.5,seed) * payload.sphere_pointer->material.diffusivity;
			current_ray.direction = current_ray.direction.reflect(new_normal);
		}
		else {
			current_colour += all_light.skycolour(current_ray.direction)*multiplier;
			break;
		}
	}
	return current_colour;
}

void render(std::list<Sphere>& objects) {
	int width = 1000;
	int height = 1000;

	//Very expensive
	int bounces = 2;
	int ray_count = 100;

	//RandomNumberGenerator rng();

	std::cout << "P3\n" << width << ' ' << height << "\n255\n";

	Light lightsource(ThreeVector(0, 8, 2), 2);
	Camera this_camera(ThreeVector(0, 4, -8), ThreeVector(0, 1, 0), ThreeVector(0, -0.5, 1));

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int rchannel = 0;
			int gchannel = 0;
			int bchannel = 0;

			double u = (x - width / 2.0) / static_cast<double>(width);
			double v = -(y - height / 2.0) / static_cast<double>(height);

			ThreeVector pixel_colour(0, 0, 0);

			if(ray_count > 1){
				for (int i = 0; i < ray_count; i++) {
					uint32_t seed = (y + i + 1) * (x + 1);
					ThreeVector ray_direction = this_camera.direction + this_camera.right * (u + (RandomDouble(seed) - 0.5) / width) + this_camera.up * (v + (RandomDouble(seed) - 0.5) / height);
					Ray ray1(this_camera.position, ray_direction);
					pixel_colour += emit_ray(ray1, objects, lightsource, bounces, seed).clamp(0, 255);
				}
				pixel_colour = pixel_colour / ray_count;
			}
			else {
				uint32_t seed = (y + 1) * (x + 1);
				ThreeVector ray_direction = this_camera.direction + this_camera.right * u + this_camera.up * v;
				Ray ray1(this_camera.position, ray_direction);
				pixel_colour = emit_ray(ray1, objects, lightsource, bounces, seed).clamp(0, 255);
			}

			rchannel = static_cast<int>(pixel_colour.x());
			gchannel = static_cast<int>(pixel_colour.y());
			bchannel = static_cast<int>(pixel_colour.z());
			std::cout << rchannel << ' ' << gchannel << ' ' << bchannel << '\n';
		}
	}
}
int main() {
	//Sphere object1(ThreeVector(2.5, 1, 0), 1, MaterialProperties(ThreeVector(200,0,0), 0.6, 0.2));
	//Sphere object2(ThreeVector(-2.5, 1, 0), 1, MaterialProperties(ThreeVector(0, 200, 0), 0.6, 0.2));
	//Sphere object3(ThreeVector(0, 1, 2), 1, MaterialProperties(ThreeVector(0,0,200), 0.6, 0.2));
	//Sphere object4(ThreeVector(0, -500, 0), 500, MaterialProperties(ThreeVector(100,100,100), 1, 0));
	//std::list<Sphere> scene = {object1,object2,object3,object4};

	Sphere object(ThreeVector(0, 0.5, 0), 3, MaterialProperties(ThreeVector(0, 0, 0), 1, 0.97));
	std::list<Sphere> scene = { object };

	auto start = std::chrono::high_resolution_clock::now();
	render(scene);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::clog << "Function took " << duration.count() << " milliseconds to execute." << std::endl;
	system("pause");
	return 0;
}
