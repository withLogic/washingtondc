/*******************************************************************************
 *
 *
 *    WashingtonDC Dreamcast Emulator
 *    Copyright (C) 2017 snickerbockers
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 ******************************************************************************/

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

#define GL3_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/gl.h>

#include "hw/pvr2/geo_buf.h"
#include "hw/pvr2/pvr2_tex_cache.h"
#include "shader.h"
#include "opengl_target.h"
#include "dreamcast.h"
#include "gfx/gfx_config.h"
#include "gfx/gfx_tex_cache.h"

#include "opengl_renderer.h"

#define POSITION_SLOT          0
#define TRANS_MAT_SLOT         1
#define BASE_COLOR_SLOT        2
#define OFFS_COLOR_SLOT        3
#define TEX_COORD_SLOT         4

static GLuint bound_tex_slot;
static GLuint tex_inst_slot;

static struct shader pvr_ta_shader;
static struct shader pvr_ta_tex_shader;

/*
 * special shader for wireframe mode that ignores vertex colors and textures
 * TODO: this shader also ignores textures.  This is not desirable since
 * textures are a separate config, but ultimately it's not that big of a deal
 * since wireframe mode always disables textures anyways
 */
static struct shader pvr_ta_no_color_shader;

static GLuint vbo, vao;

static GLuint tex_cache[PVR2_TEX_CACHE_SIZE];

static struct gfx_cfg rend_cfg;

static const GLenum tex_formats[TEX_CTRL_PIX_FMT_COUNT] = {
    [TEX_CTRL_PIX_FMT_ARGB_1555] = GL_UNSIGNED_SHORT_1_5_5_5_REV,
    [TEX_CTRL_PIX_FMT_RGB_565] = GL_UNSIGNED_SHORT_5_6_5,
    [TEX_CTRL_PIX_FMT_ARGB_4444] = GL_UNSIGNED_SHORT_4_4_4_4,
};

static const GLenum src_blend_factors[PVR2_BLEND_FACTOR_COUNT] = {
    [PVR2_BLEND_ZERO]                = GL_ZERO,
    [PVR2_BLEND_ONE]                 = GL_ONE,
    [PVR2_BLEND_OTHER]               = GL_DST_COLOR,
    [PVR2_BLEND_ONE_MINUS_OTHER]     = GL_ONE_MINUS_DST_COLOR,
    [PVR2_BLEND_SRC_ALPHA]           = GL_SRC_ALPHA,
    [PVR2_BLEND_ONE_MINUS_SRC_ALPHA] = GL_ONE_MINUS_SRC_ALPHA,
    [PVR2_BLEND_DST_ALPHA]           = GL_DST_ALPHA,
    [PVR2_BLEND_ONE_MINUS_DST_ALPHA] = GL_ONE_MINUS_DST_ALPHA
};

static const GLenum dst_blend_factors[PVR2_BLEND_FACTOR_COUNT] = {
    [PVR2_BLEND_ZERO]                = GL_ZERO,
    [PVR2_BLEND_ONE]                 = GL_ONE,
    [PVR2_BLEND_OTHER]               = GL_SRC_COLOR,
    [PVR2_BLEND_ONE_MINUS_OTHER]     = GL_ONE_MINUS_SRC_COLOR,
    [PVR2_BLEND_SRC_ALPHA]           = GL_SRC_ALPHA,
    [PVR2_BLEND_ONE_MINUS_SRC_ALPHA] = GL_ONE_MINUS_SRC_ALPHA,
    [PVR2_BLEND_DST_ALPHA]           = GL_DST_ALPHA,
    [PVR2_BLEND_ONE_MINUS_DST_ALPHA] = GL_ONE_MINUS_DST_ALPHA
};

/*
 * the PVR2 and OpenGL depth functions are inverted because PVR2's versions are
 * done based on 1 / z instead of z.  On PVR2 a closer depth-value will
 * actually be larger, and a further depth value will be smaller.  Since we
 * convert 1/z to z (in pvr2_ta.c), we also need to invert the depth comparison.
 *
 * For example, guest software which configures the depth function as
 * PVR2_DEPTH_GREATER will expect fragments with larger ("greater") depth
 * values to be in front, but after the z-component is replaced by its own
 * reciprocal, fragments with larger z-values will now have smaller z-values,
 * and fragments with smaller z-values will now have larger z-values.
 *
 * TODO: one thing I'm not sure about is whether it's coorect to convert
 * LEQUAL to GREATER, and GEQUAL to LESSER.  Mathematically these functions are
 * inversions of one another, but I'm not sure if that's what I want to do if
 * all I'm doing is accounting for the reciprocal.
 */
static const GLenum depth_funcs[PVR2_DEPTH_FUNC_COUNT] = {
    [PVR2_DEPTH_NEVER]               = GL_NEVER,
    [PVR2_DEPTH_LESS]                = GL_GEQUAL,
    [PVR2_DEPTH_EQUAL]               = GL_EQUAL,
    [PVR2_DEPTH_LEQUAL]              = GL_GREATER,
    [PVR2_DEPTH_GREATER]             = GL_LEQUAL,
    [PVR2_DEPTH_NOTEQUAL]            = GL_NOTEQUAL,
    [PVR2_DEPTH_GEQUAL]              = GL_LESS,
    [PVR2_DEPTH_ALWAYS]              = GL_ALWAYS
};

/*
 * draws the given geo_buf in whatever context is available (ie without setting
 * the shader, or the framebuffer).
 */
static void render_do_draw(struct geo_buf *geo);

// converts pixels from ARGB 4444 to RGBA 4444
static void render_conv_argb_4444(uint16_t *pixels, size_t n_pixels);

// converts pixels from ARGB 1555 to ABGR1555
static void render_conv_argb_1555(uint16_t *pixels, size_t n_pizels);

static void opengl_render_init(void);
static void opengl_render_cleanup(void);
static void opengl_renderer_update_tex(unsigned tex_obj, void const* tex_dat);
static void opengl_renderer_release_tex(unsigned tex_obj);
static void opengl_renderer_do_draw_geo_buf(struct geo_buf *geo);

struct rend_if const opengl_rend_if = {
    .init = opengl_render_init,
    .cleanup = opengl_render_cleanup,
    .update_tex = opengl_renderer_update_tex,
    .release_tex = opengl_renderer_release_tex,
    .do_draw_geo_buf = opengl_renderer_do_draw_geo_buf
};

static void opengl_render_init(void) {
    shader_load_vert_from_file(&pvr_ta_shader, "pvr2_ta_vert.glsl");
    shader_load_frag_from_file(&pvr_ta_shader, "pvr2_ta_frag.glsl");
    shader_link(&pvr_ta_shader);

    shader_load_vert_from_file_with_preamble(&pvr_ta_tex_shader,
                                             "pvr2_ta_vert.glsl",
                                             "#define TEX_ENABLE\n");
    shader_load_frag_from_file_with_preamble(&pvr_ta_tex_shader,
                                             "pvr2_ta_frag.glsl",
                                             "#define TEX_ENABLE\n");
    shader_link(&pvr_ta_tex_shader);

    shader_load_vert_from_file_with_preamble(&pvr_ta_no_color_shader,
                                             "pvr2_ta_vert.glsl",
                                             "#define COLOR_DISABLE\n");
    shader_load_frag_from_file_with_preamble(&pvr_ta_no_color_shader,
                                             "pvr2_ta_frag.glsl",
                                             "#define COLOR_DISABLE\n");
    shader_link(&pvr_ta_no_color_shader);

    bound_tex_slot = glGetUniformLocation(pvr_ta_tex_shader.shader_prog_obj,
                                          "bound_tex");
    tex_inst_slot = glGetUniformLocation(pvr_ta_tex_shader.shader_prog_obj,
                                         "tex_inst");

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenTextures(PVR2_TEX_CACHE_SIZE, tex_cache);

    unsigned tex_no;
    for (tex_no = 0; tex_no < PVR2_TEX_CACHE_SIZE; tex_no++) {
        /*
         * unconditionally set the texture wrapping mode to repeat.
         *
         * TODO: I know for sure that a lot of games need repeating texture,
         * coordinates but I don't know if there are any that need clamped
         * texture coordinates.  In the future I will need to determine if this
         * functionality exists in PVR2.
         */
        glBindTexture(GL_TEXTURE_2D, tex_no);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
}

static void opengl_render_cleanup(void) {
    glDeleteTextures(PVR2_TEX_CACHE_SIZE, tex_cache);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    shader_cleanup(&pvr_ta_no_color_shader);
    shader_cleanup(&pvr_ta_tex_shader);
    shader_cleanup(&pvr_ta_shader);

    vao = 0;
    vbo = 0;
    memset(tex_cache, 0, sizeof(tex_cache));
}

static void render_do_draw_group(struct geo_buf *geo,
                                 enum display_list_type disp_list,
                                 unsigned group_no) {
    struct poly_group *group = geo->lists[disp_list].groups + group_no;

    /*
     * TODO: currently disable color also disables textures; ideally these
     * would be two independent settings.
     */
    if (group->tex_enable && rend_cfg.tex_enable && rend_cfg.color_enable) {
        glUseProgram(pvr_ta_tex_shader.shader_prog_obj);

        if (gfx_tex_cache_get(group->tex_idx)->valid) {
            glBindTexture(GL_TEXTURE_2D, tex_cache[group->tex_idx]);
        } else {
            fprintf(stderr, "WARNING: attempt to bind invalid texture %u\n",
                    (unsigned)group->tex_idx);
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        switch (group->tex_filter) {
        case TEX_FILTER_TRILINEAR_A:
        case TEX_FILTER_TRILINEAR_B:
            printf("WARNING: trilinear filtering is not yet supported\n");
            // intentional fall-through
        case TEX_FILTER_NEAREST:
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            break;
        case TEX_FILTER_BILINEAR:
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            break;
        }

        GLenum tex_wrap_mode_gl[2];
        switch (group->tex_wrap_mode[0]) {
        case TEX_WRAP_REPEAT:
            tex_wrap_mode_gl[0] = GL_REPEAT;
            break;
        case TEX_WRAP_FLIP:
            tex_wrap_mode_gl[0] = GL_MIRRORED_REPEAT;
            break;
        case TEX_WRAP_CLAMP:
            tex_wrap_mode_gl[0] = GL_CLAMP_TO_EDGE;
            break;
        default:
            RAISE_ERROR(ERROR_INTEGRITY);
        }
        switch (group->tex_wrap_mode[1]) {
        case TEX_WRAP_REPEAT:
            tex_wrap_mode_gl[1] = GL_REPEAT;
            break;
        case TEX_WRAP_FLIP:
            tex_wrap_mode_gl[1] = GL_MIRRORED_REPEAT;
            break;
        case TEX_WRAP_CLAMP:
            tex_wrap_mode_gl[1] = GL_CLAMP_TO_EDGE;
            break;
        default:
            RAISE_ERROR(ERROR_INTEGRITY);
        }

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, tex_wrap_mode_gl[0]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, tex_wrap_mode_gl[1]);

        glUniform1i(bound_tex_slot, 0);
        glUniform1i(tex_inst_slot, group->tex_inst);
        glActiveTexture(GL_TEXTURE0);
    } else if (rend_cfg.color_enable) {
        glUseProgram(pvr_ta_shader.shader_prog_obj);
    } else {
        glUseProgram(pvr_ta_no_color_shader.shader_prog_obj);
    }

#ifdef INVARIANTS
    /*
     * this check is a little silly, but I get segfaults sometimes when
     * indexing into src_blend_factors and dst_blend_factors and I don't know
     * why.
     *
     * TODO: this was (hopefully) fixed in commit
     * 92059fe4f1714b914cec75fd2f91e676127d3097 but I am keeping the INVARIANTS
     * test here just in case.  It should be safe to delete after a couple of
     * months have gone by without this INVARIANTS test ever failing.
     */
    if (((unsigned)group->src_blend_factor >= PVR2_BLEND_FACTOR_COUNT) ||
        ((unsigned)group->dst_blend_factor >= PVR2_BLEND_FACTOR_COUNT)) {
        error_set_src_blend_factor(group->src_blend_factor);
        error_set_dst_blend_factor(group->dst_blend_factor);
        error_set_display_list_index((unsigned)disp_list);
        error_set_geo_buf_group_index(group_no);
        RAISE_ERROR(ERROR_INTEGRITY);
    }
#endif

    glBlendFunc(src_blend_factors[(unsigned)group->src_blend_factor],
                dst_blend_factors[(unsigned)group->dst_blend_factor]);


    glDepthMask(group->enable_depth_writes ? GL_TRUE : GL_FALSE);
    glDepthFunc(depth_funcs[group->depth_func]);

    /*
     * Orthographic projection.  Map all coordinates into the (-1, -1, -1) to
     * (1, 1, 1).  Anything less than -half_screen_dims or greater than
     * half_screen_dims on the x/y axes or anything not between
     * clip_min_max[0]/clip_min_max[1] on the z-axis will be clipped.  Ideally
     * nothing should be clipped on the z-axis because clip_min_max is derived
     * from the minimum and maximum depths.
     */
    GLfloat half_screen_dims[2] = {
        (GLfloat)(geo->screen_width * 0.5),
        (GLfloat)(geo->screen_height * 0.5)
    };

    /*
     * The trans_mat matrix will map z=clip_min to -1 and z=clip_max to +1.
     *
     * depending on the depth function used, this can cause fragments at the
     * extremes to not be rendered even with depth clamping enabled.  For
     * example, GL_LESS will fail anything at z=clip_min because that will get
     * mapped to -1 by this matrix, and then to +1 later in the OpenGL pipeline;
     * +1 is the furthest away value in the depth buffer, +1 will fail the depth
     * test.
     *
     * expanding clip_min and clip_max will mitigate this problem.  Obviously we
     * don't want to expand the depth range by much because that will lead to
     * precision problems, but we still want to expand it enough to make sure
     * that the minimum and maximum depth values don't get culled.
     */
    float clip_min = geo->clip_min * 1.01f;
    float clip_max = geo->clip_max * 1.01f;

    GLfloat clip_delta = clip_max - clip_min;
    GLfloat trans_mat[16] = {
        1.0 / half_screen_dims[0], 0, 0, -1,
        0, -1.0 / half_screen_dims[1], 0, 1,
        0, 0, 2.0 / clip_delta, -2.0 * clip_min / clip_delta - 1,
        0, 0, 0, 1
    };

    glUniformMatrix4fv(TRANS_MAT_SLOT, 1, GL_TRUE, trans_mat);

    // now draw the geometry itself
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(float) * group->n_verts * GEO_BUF_VERT_LEN,
                 group->verts, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(POSITION_SLOT);
    glEnableVertexAttribArray(BASE_COLOR_SLOT);
    glEnableVertexAttribArray(OFFS_COLOR_SLOT);
    glVertexAttribPointer(POSITION_SLOT, 3, GL_FLOAT, GL_FALSE,
                          GEO_BUF_VERT_LEN * sizeof(float),
                          (GLvoid*)(GEO_BUF_POS_OFFSET * sizeof(float)));
    glVertexAttribPointer(BASE_COLOR_SLOT, 4, GL_FLOAT, GL_FALSE,
                          GEO_BUF_VERT_LEN * sizeof(float),
                          (GLvoid*)(GEO_BUF_BASE_COLOR_OFFSET * sizeof(float)));
    glVertexAttribPointer(OFFS_COLOR_SLOT, 4, GL_FLOAT, GL_FALSE,
                          GEO_BUF_VERT_LEN * sizeof(float),
                          (GLvoid*)(GEO_BUF_OFFS_COLOR_OFFSET * sizeof(float)));
    if (group->tex_enable) {
        glEnableVertexAttribArray(TEX_COORD_SLOT);
        glVertexAttribPointer(TEX_COORD_SLOT, 2, GL_FLOAT, GL_FALSE,
                              GEO_BUF_VERT_LEN * sizeof(float),
                              (GLvoid*)(GEO_BUF_TEX_COORD_OFFSET * sizeof(float)));
    }
    glDrawArrays(GL_TRIANGLES, 0, group->n_verts);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glBindTexture(GL_TEXTURE_2D, 0);
}

static void render_do_draw(struct geo_buf *geo) {
}

static void opengl_renderer_update_tex(unsigned tex_obj, void const *tex_dat) {
    struct gfx_tex const *tex = gfx_tex_cache_get(tex_obj);
    GLenum format = tex->meta.pix_fmt == TEX_CTRL_PIX_FMT_RGB_565 ?
        GL_RGB : GL_RGBA;

    unsigned tex_w = 1 << tex->meta.w_shift;
    unsigned tex_h = 1 << tex->meta.h_shift;

    glBindTexture(GL_TEXTURE_2D, tex_cache[tex_obj]);
    // TODO: maybe don't always set this to 1
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    /*
     * TODO: ideally I wouldn't need to copy ARGB_4444 and ARGB_1555 into a
     * separate buffer to do the pixel conversion.  The reason I do this is that
     * the tex-dump command in the cmd thread also sees the texture data in the
     * struct gfx_tex, so I don't want to modify that.  Maybe someday I'll
     * change things to remove this mostly-unnecessary buffering...
     */
    if (tex->meta.pix_fmt == TEX_CTRL_PIX_FMT_ARGB_4444) {
        size_t n_bytes =
            sizeof(uint16_t) << (tex->meta.w_shift + tex->meta.h_shift);
        uint16_t *tex_dat_conv = (uint16_t*)malloc(n_bytes);
        if (!tex_dat_conv)
            abort();
        memcpy(tex_dat_conv, tex_dat, n_bytes);
        render_conv_argb_4444(tex_dat_conv, tex_w * tex_h);
        glTexImage2D(GL_TEXTURE_2D, 0, format, tex_w, tex_h, 0,
                     format, tex_formats[TEX_CTRL_PIX_FMT_ARGB_4444], tex_dat_conv);
        free(tex_dat_conv);
    } else if (tex->meta.pix_fmt == TEX_CTRL_PIX_FMT_ARGB_1555) {
        size_t n_bytes =
            sizeof(uint16_t) << (tex->meta.w_shift + tex->meta.h_shift);
        uint16_t *tex_dat_conv = (uint16_t*)malloc(n_bytes);
        if (!tex_dat_conv)
            abort();
        memcpy(tex_dat_conv, tex_dat, n_bytes);
        render_conv_argb_1555(tex_dat_conv, tex_w * tex_h);
        glTexImage2D(GL_TEXTURE_2D, 0, format, tex_w, tex_h, 0,
                     format, tex_formats[TEX_CTRL_PIX_FMT_ARGB_1555], tex_dat_conv);
        free(tex_dat_conv);
    } else {
        glTexImage2D(GL_TEXTURE_2D, 0, format, tex_w, tex_h, 0,
                     format, tex_formats[tex->meta.pix_fmt], tex_dat);
    }
    glBindTexture(GL_TEXTURE_2D, 0);
}

static void opengl_renderer_release_tex(unsigned tex_obj) {
    // do nothing
}

static void opengl_renderer_do_draw_geo_buf(struct geo_buf *geo) {
    gfx_config_read(&rend_cfg);

    if (!rend_cfg.wireframe) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    } else {
        glLineWidth(1);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }

    if (rend_cfg.tex_enable)
        glEnable(GL_TEXTURE_2D);
    else
        glDisable(GL_TEXTURE_2D);

    /*
     * first draw the background plane
     * TODO: I should actually draw a background plane instead
     * of just calling glClear
     */
    if (rend_cfg.bgcolor_enable) {
        glClearColor(geo->bgcolor[0], geo->bgcolor[1],
                     geo->bgcolor[2], geo->bgcolor[3]);
    } else {
        glClearColor(0.0, 0.0, 0.0, 1.0);
    }
    glDepthMask(GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (rend_cfg.depth_enable)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);

    /*
     * Strictly speaking, this isn't needed since we transform the
     * depth-component such that geo->clip_max maps to +1 and geo->clip_min
     * maps to -1, but we enable it just in case there are any floating-point
     * precision errors that push something to be greater than +1 or less
     * than -1.
     */
    glEnable(GL_DEPTH_CLAMP);

    unsigned group_no;
    enum display_list_type disp_list;
    for (disp_list = DISPLAY_LIST_FIRST; disp_list <= DISPLAY_LIST_LAST;
         disp_list++) {
        if (disp_list == DISPLAY_LIST_OPAQUE_MOD ||
            disp_list == DISPLAY_LIST_TRANS_MOD)
            continue;

        struct display_list *list = geo->lists + disp_list;

        if (rend_cfg.blend_enable) {
            if (list->blend_enable)
                glEnable(GL_BLEND);
            else
                glDisable(GL_BLEND);
        } else {
            glDisable(GL_BLEND);
        }

        if (list->n_groups && disp_list == DISPLAY_LIST_TRANS) { // TODO: only if autosort is set
            unsigned *order = malloc(sizeof(unsigned) * list->n_groups);
            rend_sort_groups(order, list->groups, list->n_groups);
            for (group_no = 0; group_no < list->n_groups; group_no++) {
                list->groups[order[group_no]].enable_depth_writes = true;
                list->groups[order[group_no]].depth_func = PVR2_DEPTH_GEQUAL;

                render_do_draw_group(geo, disp_list, order[group_no]);
            }
            free(order);
        } else {
            for (group_no = 0; group_no < list->n_groups; group_no++)
                render_do_draw_group(geo, disp_list, group_no);
        }
    }
}

static void render_conv_argb_4444(uint16_t *pixels, size_t n_pixels) {
    for (size_t pix_no = 0; pix_no < n_pixels; pix_no++, pixels++) {
        uint16_t pix_current = *pixels;
        uint16_t b = (pix_current & 0x000f) >> 0;
        uint16_t g = (pix_current & 0x00f0) >> 4;
        uint16_t r = (pix_current & 0x0f00) >> 8;
        uint16_t a = (pix_current & 0xf000) >> 12;

        *pixels = a | (b << 4) | (g << 8) | (r << 12);
    }
}

static void render_conv_argb_1555(uint16_t *pixels, size_t n_pixels) {
    for (size_t pix_no = 0; pix_no < n_pixels; pix_no++, pixels++) {
        uint16_t pix_current = *pixels;
        uint16_t b = (pix_current & 0x001f) >> 0;
        uint16_t g = (pix_current & 0x03e0) >> 5;
        uint16_t r = (pix_current & 0x7c00) >> 10;
        uint16_t a = (pix_current & 0x8000) >> 15;

        *pixels = (a << 15) | (b << 10) | (g << 5) | (r << 0);
    }
}
