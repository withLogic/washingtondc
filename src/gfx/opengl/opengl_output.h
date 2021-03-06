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

#ifndef OPENGL_OUTPUT_H_
#define OPENGL_OUTPUT_H_

/*
 * opengl_output.h: the final stage of rendering, where the framebuffer is
 * turned into and opengl texture that's rendered onto a quadrilateral
 * stretched across the screen.
 */

/*
 *
 * HOW I INTEND TO COMBINE FRAMEBUFFER RASTER GRAPHICS WITH OPENGL
 *
 * The framebuffer can reside in one of two places: in OpenGL or in WashingtonDC
 * in texture memory.  Thus, the algorithm is as follows:
 *
 * A flag will be set every time something is drawn by pvr2 via opengl; this
 * flag will indicate that the framebuffer (in PVR2 tex mem) is out of date.
 *
 * When the user reads from/writes to the framebuffer (or modifies the
 * FB_R_SOF1/2 pointer?), we will check this flag.  If it is set then
 * before the write can go through we must fetch what's in OpenGL's view of
 * the framebuffer and write it to the texture memory's framebuffer then unset
 * the flag.  This syncs OpenGL to the PVR2 framebuffer nad makes PVR2's
 * texture memory the current framebuffer.
 *
 * When there is a v-blank interrupt, we will sync OpenGL to the PVR2
 * framebuffer and then render that as a texture quad like we already do.
 * This is inefficiant and unnecessary (since games rarely/never mix PVR2 with
 * direct framebuffer access), and would be benefitted by a render-to-texture
 * approach that only syncs opengl to the PVR2 tex mem when there's a
 * read/write going on, but the lazy approach will do for now.
 *
 * If we need to do more 3d-graphics rendering while the PVR2 framebuffer is
 * the current framebuffer, then we sync the PVR2 framebuffer to OpenGL by
 * rendering it as a texture on a quad.  For this, the depth buffer will be
 * temporarily disabled, so all fragments will pass and the depth buffer will
 * not be effected.
 *
 */

#include <stdint.h>

/*
 * this gets called every time the framebuffer has a new frame to render.
 * fb_new belongs to the caller, and its contents will be copied into a new
 * storage area.
 *
 * this function is safe to call from outside of the graphics thread
 * from outside of the graphics thread, it should only be called indirectly via
 * gfx_thread_post_framebuffer.
 */
void opengl_video_new_framebuffer(uint32_t const *fb_new,
                                  unsigned fb_new_width,
                                  unsigned fb_new_height);

void opengl_video_present();

void opengl_video_output_init();
void opengl_video_output_cleanup();

#endif
