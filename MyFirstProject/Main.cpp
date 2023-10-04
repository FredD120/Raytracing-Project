#include <iostream>
#include <cmath>
#include <fstream>

class ThreeVector {
public:
	double x;
	double y;
	double z;

	ThreeVector() : x(1), y(1), z(1) {
		std::cout << "Default vector created" << std::endl;
	}

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
};

class Ray {
public:
	ThreeVector origin;
	ThreeVector direction;

	Ray(ThreeVector& o, ThreeVector& d) : origin(o), direction(d.normalise()) {}

	Ray() : origin(ThreeVector()), direction(ThreeVector().normalise()) {
		std::cout << "Default ray created" << std::endl;
	}



};

class Sphere {
public:
	ThreeVector centre;
	double radius;
	int colour;

	Sphere(ThreeVector& cen, double rad, int col) : centre(cen), radius(rad), colour(col){}

	double intersect(Ray current_ray) {
		ThreeVector origin_to_center = current_ray.origin - centre;
		double a = current_ray.direction.dot(current_ray.direction);
		double b = 2 * origin_to_center.dot(current_ray.direction);
		double c = origin_to_center.dot(origin_to_center) - radius * radius;
		double discriminant = b * b - 4 * a * c;
		if (discriminant > 0) {
			double distance1 = (-b - std::sqrt(discriminant)) / (2.0 * a);
			double distance2 = (-b + std::sqrt(discriminant)) / (2.0 * a);
			if (0 < distance1 && distance1 < distance2) {
				return distance1;
			}
			else {
				return distance2; 
			}	
		}
		else {
			return -1;
		}
	}
};


class IntArray {
public:
	static const int xsize = 10;
	static const int ysize = 10;
	int arr[xsize][ysize];
	int inital_value = 0;

	IntArray(int inval) {
		int initial_value = inval;
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				arr[i][j] = initial_value;
			}
		}
	}

	void print_array() {
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				std::cout << arr[i][j] << " ";
			}
		}
	}

	void save_array(std::string FileName) {
		std::ofstream NewFile;
		NewFile.open(FileName);
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				NewFile << arr[i][j] << " ";
			}
			NewFile << "\n" ;
		}
		NewFile.close();
	}
};

void render(Sphere object){
	std::ofstream ImageFile("C:/Users/fwdan/cImageData.txt");

	int width = 1000;
	int height = 1000;

	ThreeVector camera_position(0, 0, 0);
	ThreeVector camera_up(0, 1, 0);
	ThreeVector camera_direction(0, 0, 1);
	ThreeVector camera_right = camera_direction.cross(camera_up);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double u = (x - width / 2.0) / static_cast<double>(width);
			double v = -(y - height / 2.0) / static_cast<double>(height);
			ThreeVector ray_direction = camera_direction + camera_right.scalar_product(u) + camera_up.scalar_product(v);
			Ray ray1(camera_position, ray_direction);
			double min_distance = object.intersect(ray1);
			if (min_distance > 0) {
				int nearest_colour = object.colour;
				ImageFile << nearest_colour << " ";
			}
			else {
				ImageFile << 0 << " ";
			}
		}
		ImageFile << "\n";
	}
	ImageFile.close();
}

int main() {
	
	ThreeVector Centre1(0, 0, 3);
	double radius1 = 1;
	Sphere object1(Centre1, radius1, 1);
	render(object1);
	
	return 0;
}
