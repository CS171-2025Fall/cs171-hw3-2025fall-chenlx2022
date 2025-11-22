#include "rdr/integrator.h"

#include <omp.h>

#include "rdr/bsdf.h"
#include "rdr/camera.h"
#include "rdr/canary.h"
#include "rdr/film.h"
#include "rdr/halton.h"
#include "rdr/interaction.h"
#include "rdr/light.h"
#include "rdr/math_aliases.h"
#include "rdr/math_utils.h"
#include "rdr/platform.h"
#include "rdr/properties.h"
#include "rdr/ray.h"
#include "rdr/scene.h"
#include "rdr/sdtree.h"

RDR_NAMESPACE_BEGIN

/* ===================================================================== *
 *
 * Intersection Test Integrator's Implementation
 *
 * ===================================================================== */

void IntersectionTestIntegrator::render(ref<Camera> camera, ref<Scene> scene) {
  // Statistics
  std::atomic<int> cnt = 0;

  const Vec2i &resolution = camera->getFilm()->getResolution();
#pragma omp parallel for schedule(dynamic)
  for (int dx = 0; dx < resolution.x; dx++) {
    ++cnt;
    if (cnt % (resolution.x / 10) == 0)
      Info_("Rendering: {:.02f}%", cnt * 100.0 / resolution.x);
    Sampler sampler;
    for (int dy = 0; dy < resolution.y; dy++) {
      sampler.setPixelIndex2D(Vec2i(dx, dy));
      for (int sample = 0; sample < spp; sample++) {
        // TODO(HW3): generate #spp rays for each pixel and use Monte Carlo
        // integration to compute radiance.
        //
        // Useful Functions:
        //
        // @see Sampler::getPixelSample for getting the current pixel sample
        // as Vec2f.
        //
        // @see Camera::generateDifferentialRay for generating rays given
        // pixel sample positions as 2 floats.

        //HW3 IMPLEMENTED

        
        const Vec2f &pixel_sample = sampler.getPixelSample();
        

        auto ray = camera->generateDifferentialRay(pixel_sample.x, pixel_sample.y);


        const Vec3f &L = Li(scene, ray, sampler);
        

        assert(pixel_sample.x >= dx && pixel_sample.x <= dx + 1);
        assert(pixel_sample.y >= dy && pixel_sample.y <= dy + 1);
        camera->getFilm()->commitSample(pixel_sample, L);
      }
    }
  }
}

Vec3f IntersectionTestIntegrator::Li(
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  Vec3f color(0.0);

  // Cast a ray until we hit a non-specular surface or miss
  // Record whether we have found a diffuse surface
  bool diffuse_found = false;
  SurfaceInteraction interaction;

  for (int i = 0; i < max_depth; ++i) {
    interaction      = SurfaceInteraction();
    bool intersected = scene->intersect(ray, interaction);

    // Perform RTTI to determine the type of the surface
    bool is_ideal_diffuse =
        dynamic_cast<const IdealDiffusion *>(interaction.bsdf) != nullptr;
    bool is_perfect_refraction =
        dynamic_cast<const PerfectRefraction *>(interaction.bsdf) != nullptr;

    // Set the outgoing direction
    interaction.wo = -ray.direction;

    if (!intersected) {
      break;
    }

    if (is_perfect_refraction) {
      // We should follow the specular direction
      // TODO(HW3): call the interaction.bsdf->sample to get the new direction
      // and update the ray accordingly.
      //
      // Useful Functions:
      // @see BSDF::sample
      // @see SurfaceInteraction::spawnRay
      //
      
      //HW3 IMPLEMENTED

      interaction.bsdf->sample(interaction, sampler, nullptr);
      

      ray = interaction.spawnRay(interaction.wi);
      continue;
    }

    if (is_ideal_diffuse) {
      // We only consider diffuse surfaces for direct lighting
      diffuse_found = true;
      break;
    }

    // We simply omit any other types of surfaces
    break;
  }

  if (!diffuse_found) {
    return color;
  }

  color = directLighting(scene, interaction);
  return color;
}

Vec3f IntersectionTestIntegrator::directLighting(
    ref<Scene> scene, SurfaceInteraction &interaction) const {
  Vec3f color(0, 0, 0);
  

  //HW3 IMPLEMENTED

  if (use_area_light) {

    return directLightingAreaLight(scene, interaction);
  }
  

  Float dist_to_light = Norm(point_light_position - interaction.p);
  Vec3f light_dir     = Normalize(point_light_position - interaction.p);
  auto test_ray       = DifferentialRay(interaction.p, light_dir);

  // TODO(HW3): Test for occlusion
  //
  // You should test if there is any intersection between interaction.p and
  // point_light_position using scene->intersect. If so, return an occluded
  // color. (or Vec3f color(0, 0, 0) to be specific)
  //
  // You may find the following variables useful:
  //
  // @see bool Scene::intersect(const Ray &ray, SurfaceInteraction &interaction)
  //    This function tests whether the ray intersects with any geometry in the
  //    scene. And if so, it returns true and fills the interaction with the
  //    intersection information.
  //
  //    You can use iteraction.p to get the intersection position.
  //
  //HW3 IMPLEMENTED

  SurfaceInteraction shadow_interaction;
  

  test_ray.setTimeMax(dist_to_light - EPS);
  
  if (scene->intersect(test_ray, shadow_interaction)) {

    return Vec3f(0, 0, 0);
  }

  // Not occluded, compute the contribution using perfect diffuse diffuse model
  // Perform a quick and dirty check to determine whether the BSDF is ideal
  // diffuse by RTTI
  const BSDF *bsdf      = interaction.bsdf;
  bool is_ideal_diffuse = dynamic_cast<const IdealDiffusion *>(bsdf) != nullptr;

  if (bsdf != nullptr && is_ideal_diffuse) {
    // TODO(HW3): Compute the contribution
    //
    // You can use bsdf->evaluate(interaction) * cos_theta to approximate the
    // albedo. In this homework, we do not need to consider a
    // radiometry-accurate model, so a simple phong-shading-like model is can be
    // used to determine the value of color.

    // The angle between light direction and surface normal
    Float cos_theta =
        std::max(Dot(light_dir, interaction.normal), 0.0f);  // one-sided

    
    //HW3 IMPLEMENTED

    interaction.wi = light_dir;

    Vec3f albedo = bsdf->evaluate(interaction);
    

    Float attenuation = 0.5f / (dist_to_light * dist_to_light);
    
    color = albedo * cos_theta * point_light_flux * attenuation;
  }

  return color;
}

//HW3 IMPLEMENTED

Vec3f IntersectionTestIntegrator::directLightingAreaLight(
    ref<Scene> scene, SurfaceInteraction &interaction) const {
  Vec3f total_color(0, 0, 0);
  

  const BSDF *bsdf = interaction.bsdf;
  bool is_ideal_diffuse = dynamic_cast<const IdealDiffusion *>(bsdf) != nullptr;
  
  if (bsdf == nullptr || !is_ideal_diffuse) {
    return total_color;
  }
  

  Vec3f light_normal(0, -1, 0);  // Light facing down
  Vec3f light_tangent(1, 0, 0);  // X direction
  Vec3f light_bitangent(0, 0, 1); // Z direction
  
  int visible_samples = 0;
  

  int samples_per_dim = (int)std::sqrt((float)area_light_samples);
  int actual_samples = samples_per_dim * samples_per_dim;
  
  for (int i = 0; i < actual_samples; i++) {

    int grid_x = i % samples_per_dim;
    int grid_y = i / samples_per_dim;
    

    Float jitter_x = (Float)((i * 73) % 100) / 100.0f;
    Float jitter_y = (Float)((i * 137) % 100) / 100.0f;
    

    Float u = (grid_x + jitter_x) / samples_per_dim;
    Float v = (grid_y + jitter_y) / samples_per_dim;
    

    u = u - 0.5f;
    v = v - 0.5f;
    

    Vec3f light_sample_pos = point_light_position
        + u * area_light_size.x * light_tangent
        + v * area_light_size.y * light_bitangent;
    

    Vec3f to_light = light_sample_pos - interaction.p;
    Float dist_to_sample = Norm(to_light);
    Vec3f light_dir = to_light / dist_to_sample;  // Normalize
    

    SurfaceInteraction shadow_interaction;
    auto shadow_ray = DifferentialRay(interaction.p, light_dir);
    shadow_ray.setTimeMax(dist_to_sample - EPS);
    
    if (!scene->intersect(shadow_ray, shadow_interaction)) {

      visible_samples++;
      

      Float cos_theta = std::max(Dot(light_dir, interaction.normal), 0.0f);
      

      Float cos_theta_light = std::max(Dot(-light_dir, light_normal), 0.0f);
      

      interaction.wi = light_dir;
      Vec3f albedo = bsdf->evaluate(interaction);
      

      Float attenuation = 0.5f / (dist_to_sample * dist_to_sample);
      

      Float geometric_term = cos_theta * cos_theta_light;
      
      total_color += albedo * geometric_term * point_light_flux * attenuation;
    }
  }
  

  if (actual_samples > 0) {
    total_color = total_color / Float(actual_samples);
  }
  
  return total_color;
}

/* ===================================================================== *
 *
 * Path Integrator's Implementation
 *
 * ===================================================================== */

void PathIntegrator::render(ref<Camera> camera, ref<Scene> scene) {
  // This is left as the next assignment
  UNIMPLEMENTED;
}

Vec3f PathIntegrator::Li(
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  // This is left as the next assignment
  UNIMPLEMENTED;
}

Vec3f PathIntegrator::directLighting(
    ref<Scene> scene, SurfaceInteraction &interaction, Sampler &sampler) const {
  // This is left as the next assignment
  UNIMPLEMENTED;
}

/* ===================================================================== *
 *
 * New Integrator's Implementation
 *
 * ===================================================================== */

// Instantiate template
// clang-format off
template Vec3f
IncrementalPathIntegrator::Li<Path>(ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const;
template Vec3f
IncrementalPathIntegrator::Li<PathImmediate>(ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const;
// clang-format on

// This is exactly a way to separate dec and def
template <typename PathType>
Vec3f IncrementalPathIntegrator::Li(  // NOLINT
    ref<Scene> scene, DifferentialRay &ray, Sampler &sampler) const {
  // This is left as the next assignment
  UNIMPLEMENTED;
}

RDR_NAMESPACE_END
