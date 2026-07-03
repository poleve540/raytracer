package raytracing

import "core:fmt"
import "core:math/linalg"
import "core:math"
import "core:math/rand"

DEG2RAD :: math.RAD_PER_DEG

Vector3 :: [3]f64
Color :: distinct [3]f64

Ray :: struct
{
	origin: Vector3,
	direction: Vector3
}

Lambertian :: struct
{
	albedo: Color
}

Metal :: struct
{
	albedo: Color,
	fuzz: f64
}

Dielectric :: struct
{
	refraction_index: f64
}

Material :: union
{
	Lambertian,
	Metal,
	Dielectric
}

Sphere :: struct
{
	center: Vector3,
	radius: f64,
	mat: Material,
}

Cube :: struct
{
	center: Vector3,
	size: Vector3,
	mat: Material,
}

Cylinder :: struct
{
	base: Vector3,
	axis: Vector3,
	radius: f64,
	height: f64,
	mat: Material
}

Shape :: union
{
	Sphere,
	Cylinder
}

HitRecord :: struct
{
	point: Vector3,
	normal: Vector3,
	t: f64,
	front_face: bool,
	mat: Material
}

cylinder_hit :: proc(cylinder: Cylinder, ray: Ray, tmin, tmax: f64) -> (HitRecord, bool)
{
	a := cylinder.axis
	r := cylinder.radius
	b := cylinder.base

	oc := cylinder.base - ray.origin

	nxa := linalg.cross(ray.direction, a)

	if linalg.length(nxa) <= math.F64_EPSILON do return {}, false
	
	bnxa := linalg.dot(oc, nxa)
	discriminant := linalg.dot(nxa, nxa) * r*r - linalg.dot(a, a) * bnxa*bnxa

	if discriminant < 0 do return {}, false

	d := (linalg.dot(nxa, linalg.cross(oc, a)) - math.sqrt(discriminant)) /
	     (linalg.dot(nxa, nxa))
	t := linalg.dot(a, ray.direction*d - oc)
	if d <= tmin || d >= tmax || t < 0 || t > cylinder.height
	{
		d = (linalg.dot(nxa, linalg.cross(oc, a)) + math.sqrt(discriminant)) /
		     (linalg.dot(nxa, nxa))
		t = linalg.dot(a, ray.direction*d - oc)
		if d <= tmin || d >= tmax || t < 0 || t > cylinder.height do return {}, false
	}

	hit := HitRecord{}
	hit.point = ray.origin + ray.direction * d

	hit.normal = linalg.normalize(ray.direction*d - cylinder.axis*t - oc)

	hit.front_face = linalg.dot(ray.direction, hit.normal) < 0
	if !hit.front_face do hit.normal = -hit.normal

	hit.t = d
	hit.mat = cylinder.mat

	return hit, true
}

sphere_hit :: proc(sphere: Sphere, ray: Ray, tmin, tmax: f64) -> (HitRecord, bool)
{
	pow2 :: proc(x: $T) -> T
	{
		return linalg.pow(x, 2)
	}

	oc := sphere.center - ray.origin
	a := linalg.length2(ray.direction)
	h := linalg.dot(ray.direction, oc)
	c := linalg.length2(oc) - pow2(sphere.radius)
	discriminant := pow2(h) - a*c

	if discriminant < 0
	{
		return {}, false
	}

	sqrtd := linalg.sqrt(discriminant)
	t := (h - sqrtd) / a
	if t <= tmin || t >= tmax
	{
		t = (h + sqrtd) / a
		if t <= tmin || t >= tmax do return {}, false
	}

	hit_record := HitRecord{}
	hit_record.t = t
	hit_record.point = ray.origin + ray.direction*t
	hit_record.normal = (hit_record.point - sphere.center) / sphere.radius

	hit_record.mat = sphere.mat

	hit_record.front_face = linalg.dot(ray.direction, hit_record.normal) < 0
	if !hit_record.front_face do hit_record.normal = -hit_record.normal

	return hit_record, true
}

random_unit_vector :: proc() -> Vector3
{
	for true
	{
		p := Vector3{
			rand.float64_range(-1, 1),
			rand.float64_range(-1, 1),
			rand.float64_range(-1, 1)
		}

		len2 := linalg.length2(p)
		if len2 > math.F64_EPSILON && len2 <= 1
		{
			return p/linalg.sqrt(len2)
		}
	}

	return {}
}

reflect :: proc(v, n: Vector3) -> Vector3
{
	return v - 2*linalg.dot(v, n)*n
}

ray_color :: proc(ray: Ray, shapes: []Shape, depth: int) -> Color
{
	if depth <= 0 do return Color{}

	hit := HitRecord{}
	hit_anything := false
	closest_so_far := math.inf_f64(1)
	for shape in shapes
	{
		switch s in shape
		{
		case Sphere:
			hit = sphere_hit(s, ray, 0.001, closest_so_far) or_continue
		case Cylinder:
			hit = cylinder_hit(s, ray, 0.001, closest_so_far) or_continue
		}

		hit_anything = true
		closest_so_far = hit.t
	}

	if hit_anything
	{
		switch mat in hit.mat
		{
		case Dielectric:
			direction := Vector3{}

			refraction_index := hit.front_face ? 1/mat.refraction_index : mat.refraction_index
			unit_dir := linalg.normalize(ray.direction)
			cos_theta := min(linalg.dot(-unit_dir, hit.normal), 1.0)
			sin_theta := math.sqrt(1.0 - cos_theta*cos_theta)

			r0 := (1 - refraction_index) / (1 + refraction_index)
			r0 = r0 * r0
			reflectance := r0 + (1-r0) * math.pow((1 - cos_theta), 5)

			if refraction_index * sin_theta > 1 || reflectance > rand.float64()
			{
				direction = reflect(unit_dir, hit.normal)
			}
			else
			{
				r_out_perp := refraction_index * (unit_dir + cos_theta*hit.normal)
				r_out_parallel := -math.sqrt(math.abs(1.0 - linalg.length2(r_out_perp))) * hit.normal
				direction = r_out_perp + r_out_parallel
			}

			return Color(1) * ray_color({hit.point, direction}, shapes, depth-1)

		case Lambertian:
			reflected_direction := random_unit_vector() + hit.normal
			if linalg.all(linalg.less_than(reflected_direction, Vector3(math.F64_EPSILON)))
			{
				reflected_direction = hit.normal
			}
			return mat.albedo * ray_color({hit.point, reflected_direction}, shapes, depth-1)

		case Metal:
			reflected := reflect(ray.direction, hit.normal)
			reflected = linalg.normalize(reflected) + (mat.fuzz * random_unit_vector())
			if linalg.dot(reflected, hit.normal) > 0
			{
				return mat.albedo * ray_color({hit.point, reflected}, shapes, depth-1)
			}
			else
			{
				return Color{}
			}
		}
	}

	a := (linalg.normalize(ray.direction).y + 1)/2
	return linalg.lerp(Color{1, 1, 1}, Color{0.5, 0.7, 1}, a)
}

main :: proc()
{
	IMAGE_WIDTH :: 400
	IMAGE_HEIGHT :: 225

	fov: f64 = 40*DEG2RAD

	lookfrom := Vector3{-2, 2, 1}
	lookat := Vector3{0, 0, -1}
	vup := Vector3{0, 1, 0}

	camera_center := lookfrom

	defocus_angle: f64 = 0
	focus_dist: f64 = 10

	w := linalg.normalize(lookfrom - lookat)
	u := linalg.normalize(linalg.cross(vup, w))
	v := linalg.cross(w, u)

	viewport_height: f64 = 2 * math.tan(fov/2) * focus_dist
	viewport_width := viewport_height * f64(IMAGE_WIDTH)/IMAGE_HEIGHT
	viewport_u := viewport_width * u
	viewport_v := viewport_height * -v

	pixel_delta_u := viewport_u / IMAGE_WIDTH
	pixel_delta_v := viewport_v / IMAGE_HEIGHT

	viewport_upperleft := camera_center - (focus_dist * w) - viewport_u/2 - viewport_v/2
	pixel00_loc := viewport_upperleft + 0.5 * (pixel_delta_u + pixel_delta_v)

	defocus_radius := focus_dist * math.tan(math.to_radians(defocus_angle / 2))
	defocus_disk_u := u * defocus_radius
	defocus_disk_v := v * defocus_radius

	fmt.println("P3")
	fmt.println(IMAGE_WIDTH, IMAGE_HEIGHT)
	fmt.println("255")

	// --- 3 spheres ---
	// spheres := []Shape{
	// 	Sphere{{0, -100.5, -1}, 100, Lambertian{{0.8, 0.8, 0}}},
	// 	Sphere{{0, 0, -1.2}, 0.5, Lambertian{{0.1, 0.2, 0.5}}},
	// 	Sphere{{-1, 0, -1}, 0.5, Dielectric{1.5}},
	// 	Sphere{{-1, 0, -1}, 0.4, Dielectric{1/1.5}},
	// 	Sphere{{1, 0, -1}, 0.5, Metal{{0.8, 0.6, 0.2}, 1}},
	// }

	// shapes := []Shape{
	// 	Sphere{{0, -100.5, -1}, 100, Lambertian{{0.8, 0.8, 0}}},
	// 	Sphere{{0, 0, -1.2}, 0.5, Lambertian{{0.1, 0.2, 0.5}}},
	// }
	//
	// cylinder := Cylinder{{0, -0.5, -1.2}, {0, 1, 0}, 0.5, 1.0, Lambertian{{0.1, 0.2, 0.5}}}
	// shapes[1] = cylinder

	//--- Weiner ---
	// spheres := make([dynamic]Shape)
	// append(&spheres, Sphere{{0, -1000, 0}, 1000, Lambertian{{0.8, 0.5, 0.5}}})
	//
	// skin := Lambertian{{232.0/256, 174.0/256, 120.0/256}}
	// append(&spheres, Sphere{{0, 0.2, 0}, 0.2, skin})
	// append(&spheres, Sphere{{-0.3, 0.2, 0}, 0.2, skin})
	// append(&spheres, Sphere{{-0.6, 0.2, 0}, 0.2, skin})
	//
	// // tip := Lambertian{{239.0/256, 160.0/256, 230.0/256}}
	// tipin := Dielectric{1.0/1.5}
	// tipout := Dielectric{1.5}
	// append(&spheres, Sphere{{-0.9, 0.21, 0}, 0.15, tipin})
	// append(&spheres, Sphere{{-0.9, 0.21, 0}, 0.2, tipout})
	//
	// ball := Metal{0.8, 0}
	// append(&spheres, Sphere{{0.4, 0.2, 0.4}, 0.4, ball})
	// append(&spheres, Sphere{{0.4, 0.2, -0.4}, 0.4, ball})

	// --- Cover ---
	spheres: [dynamic]Shape
	for a in -11..<11
	{
		for b in -11..<11
		{
			choose_mat := rand.float64()
			center := Vector3{cast(f64)a + 0.9*rand.float64(), 0.2, cast(f64)b + 0.9*rand.float64()}
			if linalg.distance(Vector3{4, 0.2, 0}, center) > 0.9
			{
				switch
				{
				case choose_mat < 0.8:
					albedo := Color{rand.float64(),rand.float64(),rand.float64()}
					append(&spheres, Sphere{center, 0.2, Lambertian{albedo * albedo}})

				case choose_mat < 0.95:
					albedo := Color{rand.float64_range(0.5,1),rand.float64_range(0.5,1),rand.float64_range(0.5,1)}
					fuzz := rand.float64_range(0, 0.5)
					append(&spheres, Sphere{center, 0.2, Metal{albedo, fuzz}})

				case:
					append(&spheres, Sphere{center, 0.2, Dielectric{1.5}})
				}
			}
		}
	}
	append(&spheres, Sphere{{0, 1, 0}, 1, Dielectric{1.5}})
	append(&spheres, Sphere{{-4, 1, 0}, 1, Lambertian{{0.4, 0.2, 0.1}}})
	append(&spheres, Sphere{{4, 1, 0}, 1, Metal{{0.7, 0.6, 0.5}, 0}})

	NUM_SAMPLES_PER_PIXEL :: 50
	MAX_DEPTH :: 50

	for y in 0..<IMAGE_HEIGHT
	{
		fmt.eprintfln("Scanlines remaining: %v", IMAGE_HEIGHT-y)

		for x in 0..<IMAGE_WIDTH
		{
			color := Color{}
			for sample in 0..<NUM_SAMPLES_PER_PIXEL
			{
				offset := Vector3{rand.float64()-0.5, rand.float64()-0.5, 0}
				pixel_sample := pixel00_loc + ((cast(f64)x + offset.x) * pixel_delta_u) + ((cast(f64)y + offset.y) * pixel_delta_v)

				p := Vector3{}
				for true
				{
					p = Vector3{rand.float64_range(-1, 1), rand.float64_range(-1, 1), 0}
					len2 := linalg.length2(p)
					if len2 < 1 do break
				}

				ray := Ray{}
				ray.origin = defocus_angle > 0 \
						 ? camera_center + (p.x * defocus_disk_u) + (p.y * defocus_disk_v) \
						 : camera_center
				ray.direction = pixel_sample-ray.origin

				color += ray_color(ray, spheres[:], MAX_DEPTH)
			}

			color /= NUM_SAMPLES_PER_PIXEL

			for &c in color
			{
				if c > 0 do c = math.sqrt(c)
			}

			color = linalg.clamp(color, 0, 0.999)
			fmt.println(int(256 * color.r), int(256 * color.g), int(256 * color.b))
		}
	}
}
