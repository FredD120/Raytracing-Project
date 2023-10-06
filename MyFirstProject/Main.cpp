#include <iostream>
#include <cmath>
#include <fstream>
#include <list>
#include <limits>


class ThreeVector {
public:
	double x;
	double y;
	double z;

	ThreeVector() : x(1), y(1), z(1) {}

	ThreeVector(double entry1, double entry2, double entry3) :
		x(entry1), y(entry2), z(entry3) {}


	void print_vector() {
		std::cout << x << "\n" << y << "\n" << z << std::endl;
	}
	
	double dot(ThreeVector other) {
		return x * other.x + y * other.y + z * other.z;
	}

	ThreeVector cross(ThreeVector other) {
		return ThreeVector(y * other.z - z * other.y, z* other.x - x * other.z, x* other.y - y * other.x);
	}

	ThreeVector normalise() {
		double mag = magnitude();
		if (mag == 0){
			std::cout << "Magnitude zero vector cannot be normalised" << std::endl;
			return *this;
		} else {
			return ThreeVector(x / mag, y / mag, z / mag);
		}
	}

	ThreeVector scalar_product(double alpha) {
		return ThreeVector(x * alpha, y * alpha, z * alpha);
	}

	double magnitude() {
		return std::sqrt(x*x + y*y + z*z);
	}

	ThreeVector operator-(ThreeVector other) {
		return ThreeVector(x - other.x, y - other.y, z - other.z);
	}

	ThreeVector operator+(ThreeVector other) {
		return ThreeVector(x + other.x, y + other.y, z + other.z);
	}
	ThreeVector invert() {
		return ThreeVector(-x, -y, -z);
	}
};

class Ray {
public:
	ThreeVector origin;
	ThreeVector direction;

	Ray(ThreeVector& o, ThreeVector& d) : origin(o), direction(d) {}

	Ray() : origin(ThreeVector()), direction(ThreeVector()) {}
};

class Sphere {
public:
	ThreeVector centre;
	double radius;
	ThreeVector colour;

	Sphere(ThreeVector& cen, double rad, ThreeVector col) : centre(cen), radius(rad), colour(col){}

	Sphere() : centre(ThreeVector()), radius(1), colour(ThreeVector()) {}

	double intersect(Ray current_ray) {
		ThreeVector origin_to_center = current_ray.origin - centre;
		double a = current_ray.direction.dot(current_ray.direction);
		double b = 2 * origin_to_center.dot(current_ray.direction);
		double c = origin_to_center.dot(origin_to_center) - radius * radius;
		double discriminant = b * b - 4 * a * c;
		if (discriminant > 0) {
			double distance1 = (-b - std::sqrt(discriminant)) / (2.0 * a);
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

	ThreeVector normal(ThreeVector position) {
		return (position - centre).normalise();
	}
};

void render(std::list<Sphere> objects){
	std::ofstream ImageFileR("C:/Users/fwdan/C++ Raytracing/ImageDataR.txt");
	std::ofstream ImageFileG("C:/Users/fwdan/C++ Raytracing/ImageDataG.txt");
	std::ofstream ImageFileB("C:/Users/fwdan/C++ Raytracing/ImageDataB.txt");

	int width = 1000;
	int height = 1000;

	double lightsource_brightness = 0.8;
	double ambient_brightness = 1-lightsource_brightness;
	ThreeVector lightsource_position(1.5,1,-2);

	ThreeVector camera_position(0, 0, -5);
	ThreeVector camera_up(0, 1, 0);
	ThreeVector camera_direction(0, -0.2, 1);
	ThreeVector camera_right = camera_direction.cross(camera_up);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double u = (x - width / 2.0) / static_cast<double>(width);
			double v = -(y - height / 2.0) / static_cast<double>(height);
			ThreeVector ray_direction = camera_direction + camera_right.scalar_product(u) + camera_up.scalar_product(v);
			Ray ray1(camera_position, ray_direction);

			double min_distance = std::numeric_limits<double>::max();
			Sphere* current_object_pointer = nullptr;

			for (Sphere& object : objects) {
				double current_distance = object.intersect(ray1);
				if (current_distance < min_distance) {
					min_distance = current_distance;
					current_object_pointer = &object;
				}
			}			
			if (current_object_pointer != nullptr) {
				ThreeVector intersect_pos = camera_position + ray_direction.scalar_product(min_distance);
				ThreeVector light_to_sphere = (lightsource_position - intersect_pos).normalise();
				double lightsource_reflection = lightsource_brightness * current_object_pointer->normal(intersect_pos).dot(light_to_sphere);
				ThreeVector nearest_colour = ThreeVector();
				if (lightsource_reflection>0){
					nearest_colour = current_object_pointer->colour.scalar_product(ambient_brightness + lightsource_reflection);
				}
				else {
					nearest_colour = current_object_pointer->colour.scalar_product(ambient_brightness);
				}
				ImageFileR << static_cast<int>(nearest_colour.x) << " ";
				ImageFileG << static_cast<int>(nearest_colour.y) << " ";
				ImageFileB << static_cast<int>(nearest_colour.z) << " ";
			}
			else {
				ImageFileR << 0 << " ";
				ImageFileG << 0 << " ";
				ImageFileB << 0 << " ";
			}
		}
		ImageFileR << "\n";
		ImageFileG << "\n";
		ImageFileB << "\n";
	}
	ImageFileR.close();
	ImageFileG.close();
	ImageFileB.close();
}

int main() {
	
	ThreeVector Centre1(1, 0.5, 6);
	double radius1 = 1;
	ThreeVector colour1(200, 0, 200);

	ThreeVector Centre2(-1, -0.5, 4);
	double radius2 = 1.4;
	ThreeVector colour2(0, 150, 250);

	ThreeVector Centre3(0, -100, 3);
	double radius3 = 95;
	ThreeVector colour3(50, 250, 50);

	ThreeVector light(2, -1, 3);
	double light_radius = 0.5;
	ThreeVector light_colour(250, 250, 250);

	Sphere object1(Centre1, radius1, colour1);
	Sphere object2(Centre2, radius2, colour2);
	Sphere object3(Centre3, radius3, colour3);
	Sphere lightsource(light, light_radius, light_colour);
	std::list<Sphere> scene = {object1,object2,object3,lightsource};
	render(scene);
	
	return 0;
}
