#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

const double PI = 3.14159265358979323846;
const double INF = 1e20;
const double EPS = 1e-6;

// *** その他の関数 ***
inline double rand01() { return (double)rand()/RAND_MAX; }

// *** データ構造 ***
struct Vec {
	double x, y, z;
	Vec(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
	inline Vec operator+(const Vec &b) const {return Vec(x + b.x, y + b.y, z + b.z);}
	inline Vec operator-(const Vec &b) const {return Vec(x - b.x, y - b.y, z - b.z);}
	inline Vec operator*(const double b) const {return Vec(x * b, y * b, z * b);}
	inline Vec operator/(const double b) const {return Vec(x / b, y / b, z / b);}
	inline const double LengthSquared() const { return x*x + y*y + z*z; }
	inline const double Length() const { return sqrt(LengthSquared()); }
};
inline Vec operator*(double f, const Vec &v) { return v * f; }
inline Vec Normalize(const Vec &v) { return v / v.Length(); }
// 要素ごとの積をとる
inline const Vec Multiply(const Vec &v1, const Vec &v2) {
	return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline const double Dot(const Vec &v1, const Vec &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline const Vec Cross(const Vec &v1, const Vec &v2) {
	return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
}
typedef Vec Color;
const Color BackgroundColor(0.0, 0.0, 0.0);

struct Ray {
	Vec org, dir;
	Ray(const Vec org_, const Vec &dir_) : org(org_), dir(dir_) {}
};

enum ReflectionType {
	DIFFUSE,    // 完全拡散面。いわゆるLambertian面。
};

struct Rectangle {
	Color emission, color;
	Vec a, b, p0;
	Vec normal;
	
	std::vector<Color> patch; // パッチのラジオシティ
	int a_num, b_num;
	double a_len, b_len;

	inline Color sample(int ia, int ib) const {
		if (ia < 0) ia = 0;
		if (ia >= a_num) ia = a_num - 1;
		if (ib < 0) ib = 0;
		if (ib >= b_num) ib = b_num - 1;
		return patch[ia * b_num + ib];
	}

	void divide_patchs(const int a_num_, const int b_num_) {
		a_num = a_num_;
		b_num = b_num_;
		patch.clear();
		patch.resize(a_num * b_num);
	}

	Rectangle(const Vec p0_, const Vec &a_, const Vec &b_, const Color &emission_, const Color &color_) :
	a(a_), b(b_), p0(p0_), emission(emission_), color(color_) {
		normal = Cross(a, b);
		normal = Normalize(normal);

		a_len = a.Length();
		b_len = b.Length();
	}
	inline const double intersect(const Ray &ray) {
		const double t = Dot(p0 - ray.org, normal) / Dot(ray.dir, normal);
		if (t <= EPS)
			return 0.0;

		Vec p = ray.org + t * ray.dir;
		Vec d = p - p0;
		const double ddota = Dot(d, a);
		if (ddota < 0.0 || ddota > a.LengthSquared())
			return 0.0;

		const double ddotb = Dot(d, b);
		if (ddotb < 0.0 || ddotb > b.LengthSquared())
			return 0.0;

		return t;
	}

};

// *** レンダリングするシーンデータ ****
// from small ppt
Rectangle recs[] = {
	
	Rectangle(Vec(0.0, 0.0, 0.0),     Vec(100.0, 0.0, 0.0),  Vec(0.0, 80.0, 0.0),   Color(), Color(0.75, 0.75, 0.75)), // 奥
	Rectangle(Vec(0.0, 0.0, 170.0),   Vec(100.0, 0.0, 0.0),  Vec(0.0, 0.0, -170.0), Color(), Color(0.75, 0.75, 0.75)), // 床
	Rectangle(Vec(0.0, 80.0, 0.0),    Vec(100.0, 0.0, 0.0),  Vec(0.0, 0.0, 170.0),  Color(), Color(0.75, 0.75, 0.75)), // 天井
	Rectangle(Vec(0.0, 0.0, 170.0),   Vec(0.0, 0.0, -170.0), Vec(0.0, 80.0, 0.0),   Color(), Color(0.75, 0.25, 0.25)), // 左
	Rectangle(Vec(100.0, 0.0, 0.0),   Vec(0.0, 0.0, 170.0),  Vec(0.0, 80.0, 0.0),   Color(), Color(0.25, 0.25, 0.75)), // 右
	Rectangle(Vec(100.0, 0.0, 170.0), Vec(-100.0, 0.0, 0.0), Vec(0.0, -80.0, 0.0),  Color(), Color()), // 手前
	
	Rectangle(Vec(40.0, 79.99, 65.0), Vec(20.0, 0.0, 0.0), Vec(0.0, 0.0, 20.0), Color(6,6,6) / PI, Color(0.75, 0.75, 0.75)), // 照明

	Rectangle(Vec(30.0, 0.0, 100.0),   Vec(0.0, 0.0, -20.0), Vec(0.0, 40.0, 0.0),   Color(), Color(0.75, 0.75, 0.75)), // 箱右
	Rectangle(Vec(10.0, 0.0, 80.0),   Vec(0.0, 0.0, 20.0),  Vec(0.0, 40.0, 0.0),   Color(), Color(0.75, 0.75, 0.75)), // 箱左
	Rectangle(Vec(10.0, 0.0, 100.0),     Vec(20.0, 0.0, 0.0),  Vec(0.0, 40.0, 0.0),   Color(), Color(0.75, 0.75, 0.75)), // 箱手前
	Rectangle(Vec(30.0, 0.0, 80.0), Vec(-20.0, 0.0, 0.0), Vec(0.0, -40.0, 0.0),  Color(), Color(0.75, 0.75, 0.75)), // 箱奥
	Rectangle(Vec(10.0, 40.0, 100.0),   Vec(20.0, 0.0, 0.0),  Vec(0.0, 0.0, -20.0), Color(), Color(0.75, 0.75, 0.75)), // 箱天井
};


// *** レンダリング用関数 ***
// シーンとの交差判定関数
inline bool intersect_scene(const Ray &ray, double *t, int *id, Vec *normal) {
	const double n = sizeof(recs) / sizeof(Rectangle);
	*t  = INF;
	*id = -1;
	for (int i = 0; i < int(n); i ++) {
		double d = recs[i].intersect(ray);
		if (d > 0.0 && d < *t) {
			*t  = d;
			*id = i;
			*normal = recs[i].normal;
		}
	}
	return *t < INF;
}

static double *form_factor;
static int patch_num = 0;

void calc_form_factor(const int a_div_num, const int b_div_num, const int mc_sample) {
	const double n = sizeof(recs) / sizeof(Rectangle);
	for (int i = 0; i < int(n); i ++) {
		recs[i].divide_patchs(a_div_num, b_div_num); // とりあえず適当に分割
		patch_num += recs[i].a_num * recs[i].b_num;
	}

	std::cout << "patch num: " << patch_num << " form factor num:" << patch_num * patch_num << std::endl;

	form_factor = new double[patch_num * patch_num];
	memset(form_factor, 0.0, sizeof(double) * patch_num * patch_num);
	
//#pragma omp parallel for schedule(dynamic, 1) num_threads(3)
	for (int i = 0; i < int(n); i ++) {
		srand(i * i * i + i);
		int patch_i = 0;
		for (int k = 0; k < i; k ++)
			patch_i += recs[k].a_num * recs[k].b_num;
		std::cout << i << std::endl;
		for (int ia = 0; ia < recs[i].a_num; ia ++) {
			std::cout << "*";
			for (int ib = 0; ib < recs[i].b_num; ib ++) {
				// 面(i)上の、パッチ(ia, ib)
				// form_factor[patch_i * patch_num + patch_j] = F_ij
				const Vec normal_i = recs[i].normal;
				const double Ai = Cross(recs[i].a / recs[i].a_num, recs[i].b / recs[i].b_num).Length(); // パッチ(ia, ib)の面積

				// 相手の面
				int patch_j = 0;
				for (int j = 0; j < int(n); j ++) {
					const Vec normal_j = recs[j].normal;
					for (int ja = 0; ja < recs[j].a_num; ja ++) {
						for (int jb = 0; jb < recs[j].b_num; jb ++) {
							const double Aj = Cross(recs[j].a / recs[j].a_num, recs[j].b / recs[j].b_num).Length(); // パッチ(ja, jb)の面積
							if (i != j) {
								// フォームファクター計算
								// モンテカルロ積分
								// 等間隔にサンプリング
								double F = 0;
								const int Ni = mc_sample, Nj = mc_sample;
								for (int ias = 0; ias < Ni; ias ++) {
									for (int ibs = 0; ibs < Ni; ibs ++) {
										for (int jas = 0; jas < Nj; jas ++) {
											for (int jbs = 0; jbs < Nj; jbs ++) {
												const double u0 = (double)(ias + 0.5) / Ni, u1 = (double)(ibs + 0.5) / Ni;
												const double u2 = (double)(jas + 0.5) / Nj, u3 = (double)(jbs + 0.5) / Nj;
												const Vec xi = recs[i].p0 + recs[i].a * ((double)(ia + u0) / recs[i].a_num) + recs[i].b * ((double)(ib + u1) / recs[i].b_num);
												const Vec xj = recs[j].p0 + recs[j].a * ((double)(ja + u2) / recs[j].a_num) + recs[j].b * ((double)(jb + u3) / recs[j].b_num);

												// V(x, y)
												const Vec ij = Normalize(xj - xi);
												double t; // レイからシーンの交差位置までの距離
												int id;   // 交差したシーン内オブジェクトのID
												Vec normal; // 交差位置の法線
												if (intersect_scene(Ray(xi, ij), &t, &id, &normal) && id != j) {
													continue;
												}

												const double d0 = Dot(normal_i, ij);
												const double d1 = Dot(normal_j, -1.0 * ij);

												if (d0 > 0.0 && d1 > 0.0) {
													const double K = d0 * d1 / (PI * (xj - xi).LengthSquared());
													const double pdf = (1.0 / Ai) * (1.0 / Aj);
													F += K / pdf;
												}

											}
										}
									}
								} 
								F /= (Ni + 1) * (Ni + 1) * (Nj + 1) * (Nj + 1) * Ai;
								if (F > 1.0)
									F = 1.0;
								form_factor[patch_i * patch_num + patch_j] = F;
							}
							patch_j ++;
						}
					}
				}
				patch_i ++;
			}
		}
	}
}

void calc_radiosity(const int iteration) {
	// ガウス・ザイデル法で連立一次方程式（ラジオシティ方程式）を解く
	const double n = sizeof(recs) / sizeof(Rectangle);
	int patch_i = 0;
	for (int i = 0; i < int(n); i ++) {
		for (int ia = 0; ia < recs[i].a_num; ia ++) {
			for (int ib = 0; ib < recs[i].b_num; ib ++) {
				// 面(i)上の、パッチ(ia, ib)
				Color B;

				// 相手の面
				int patch_j = 0;
				for (int j = 0; j < int(n); j ++) {
					for (int ja = 0; ja < recs[j].a_num; ja ++) {
						for (int jb = 0; jb < recs[j].b_num; jb ++) {
							const double Fij = form_factor[patch_i * patch_num + patch_j];
							if (Fij > 0.0)
								B = B + Fij * recs[j].patch[ja * recs[j].b_num + jb];
							patch_j ++;
						}
					}
				}
				B = Multiply(recs[i].color, B) + recs[i].emission;

				recs[i].patch[ia * recs[i].b_num + ib] = B;
				patch_i ++;
			}
		}
	}
}

// http://www.paulinternet.nl/?page=bicubic より拝借
Color cubicInterpolate (Color p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

Color bicubicInterpolate (Color p[4][4], double x, double y) {
	Color arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}

// ray方向からの放射輝度を求める
Color radiance(const Ray &ray, const int depth, bool interpolation = true) {
	double t; // レイからシーンの交差位置までの距離
	int id;   // 交差したシーン内オブジェクトのID
	Vec normal; // 交差位置の法線
	if (!intersect_scene(ray, &t, &id, &normal))
		return BackgroundColor;

	const Rectangle &obj = recs[id];
	const Vec hitpoint = ray.org + t * ray.dir; // 交差位置

	// シーンを構成するジオメトリは全て完全拡散面と仮定
	switch (DIFFUSE) {
	case DIFFUSE: {
		const Vec v = hitpoint - obj.p0;
		const double a_len = Dot(v, Normalize(obj.a));
		const double b_len = Dot(v, Normalize(obj.b));

		double da = obj.a_num * a_len / obj.a_len;
		double db = obj.b_num * b_len / obj.b_len;

		int ia = int(da); if (ia >= obj.a_num) ia --;
		int ib = int(db); if (ib >= obj.b_num) ib --;

		// バイキュービック補間
		if (interpolation) {
			Color c[4][4];

			int ia = int(da - 0.5);
			int ib = int(db - 0.5);

			for (int i = 0; i < 4; i ++) {
				for (int j = 0; j < 4; j ++) {
					c[i][j] = obj.sample(ia + i - 1, ib + j - 1);
				}
			}
			
			int ia0 = int(da - 0.5);
			int ib0 = int(db - 0.5); 
			double dx = da - ia0 - 0.5;
			double dy = db - ib0 - 0.5;
			return PI * bicubicInterpolate(c, dx, dy);
		}
		// 完全拡散面なのでPI倍
		return PI * obj.patch[ia * obj.b_num + ib];
	} break;
	}
}


// *** .hdrフォーマットで出力するための関数 ***
struct HDRPixel {
	unsigned char r, g, b, e;
	HDRPixel(const unsigned char r_ = 0, const unsigned char g_ = 0, const unsigned char b_ = 0, const unsigned char e_ = 0) :
	r(r_), g(g_), b(b_), e(e_) {};
	unsigned char get(int idx) {
		switch (idx) {
		case 0: return r;
		case 1: return g;
		case 2: return b;
		case 3: return e;
		} return 0;
	}

};

// doubleのRGB要素を.hdrフォーマット用に変換
HDRPixel get_hdr_pixel(const Color &color) {
	double d = std::max(color.x, std::max(color.y, color.z));
	if (d <= 1e-32)
		return HDRPixel();
	int e;
	double m = frexp(d, &e); // d = m * 2^e
	d = m * 256.0 / d;
	return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
}

// 書き出し用関数
void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cerr << "Error: " << filename << std::endl;
		return;
	}
	// .hdrフォーマットに従ってデータを書きだす
	// ヘッダ
	unsigned char ret = 0x0a;
	fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
	fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
	fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
	fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

	// 輝度値書き出し
	fprintf(fp, "-Y %d +X %d%c", height, width, ret);
	for (int i = height - 1; i >= 0; i --) {
		std::vector<HDRPixel> line;
		for (int j = 0; j < width; j ++) {
			HDRPixel p = get_hdr_pixel(image[j + i * width]);
			line.push_back(p);
		}
		fprintf(fp, "%c%c", 0x02, 0x02);
		fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
		for (int i = 0; i < 4; i ++) {
			for (int cursor = 0; cursor < width;) {
				const int cursor_move = std::min(127, width - cursor);
				fprintf(fp, "%c", cursor_move);
				for (int j = cursor;  j < cursor + cursor_move; j ++) {
					fprintf(fp, "%c", line[j].get(i));
				}
				cursor += cursor_move;
			}
		}
	}

	fclose(fp);
}

int main(int argc, char **argv) {
	int width = 640;
	int height = 480;
	int samples = 16;

	// カメラ位置
	Ray camera(Vec(50.0, 52.0, 295.6), Normalize(Vec(0.0, -0.042612, -1.0)));
	// シーン内でのスクリーンのx,y方向のベクトル
	Vec cx = Vec(width * 0.5135 / height);
	Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;
	Color *image = new Color[width * height];
	Color *image_interpolated = new Color[width * height];

	// フォームファクター計算
	// 16x16に分割、フォームファクターの計算に(3x3)^2サンプル使ってモンテカルロ積分する
	calc_form_factor(16, 16, 3);

	// ガウス=ザイデル法でラジオシティ方程式解く
	// 32反復
	for (int i = 0; i < 32; i ++) {
		std::cout << i << " ";
		calc_radiosity(i);
	}

//#pragma omp parallel for schedule(dynamic, 1) num_threads(3)
	for (int y = 0; y < height; y ++) {
		std::cerr << "Rendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%" << std::endl;
		srand(y * y * y);
		for (int x = 0; x < width; x ++) {
			int image_index = y * width + x;	
			image[image_index] = Color();
			image_interpolated[image_index] = Color();

			// 2x2のサブピクセルサンプリング
			for (int sy = 0; sy < 2; sy ++) {
				for (int sx = 0; sx < 2; sx ++) {
					Color accumulated_radiance = Color();
					Color accumulated_radiance2 = Color();
					// 一つのサブピクセルあたりsamples回サンプリングする
					for (int s = 0; s < samples; s ++) {
						// テントフィルターによってサンプリング
						// ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
						const double r1 = 2.0 * rand01(), dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
						const double r2 = 2.0 * rand01(), dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
						Vec dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
								  cy * (((sy + 0.5 + dy) / 2.0 + y) / height- 0.5) + camera.dir;
						accumulated_radiance = accumulated_radiance + 
							radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 0, false) / samples;
						accumulated_radiance2 = accumulated_radiance2 + 
							radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 0, true) / samples;
					}
					image[image_index] = image[image_index] + accumulated_radiance;
					image_interpolated[image_index] = image_interpolated[image_index] + accumulated_radiance2;

				}
			}
		}
	}
	
	// .hdrフォーマットで出力
	save_hdr_file(std::string("image.hdr"), image, width, height);
	save_hdr_file(std::string("image_interpolated.hdr"), image_interpolated, width, height);
}
