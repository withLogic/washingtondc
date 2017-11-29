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

/*
 * common structures and interfaces that are used between the rendering code
 * and the non-rendering code.  Everything here pertains exclusively to the
 * gfx_thread.
 */

#ifndef REND_COMMON_H_
#define REND_COMMON_H_

#include "hw/pvr2/geo_buf.h"

struct rend_if {
    void (*init)(void);

    void (*cleanup)(void);

    /*
     * called to notify the renderer that it needs to update the given
     * texture from the texture cache
     */
    void (*update_tex)(unsigned tex_obj, void const*);

    /*
     * called to notify the renderer that it needs to release the resources
     * associated with the given texture.
     */
    void (*release_tex)(unsigned tex_obj);

    /*
     * render the contents of the given geo_buf
     *
     * TODO: split this up into multiple functions (apply config, apply
     * graphics state, render list, etc).
     */
    void (*do_draw_geo_buf)(struct geo_buf *geo);
};

// initialize and clean up the graphics renderer
void rend_init(void);
void rend_cleanup(void);

// tell the renderer to update the given texture from the cache
void rend_update_tex(unsigned tex_no, void const *tex_dat);

// tell the renderer to release the given texture from the cache
void rend_release_tex(unsigned tex_no);

// draw the given geo_buf
void rend_do_draw_geo_buf(struct geo_buf *geo);

// this should only be called from the gfx_thread
void rend_draw_geo_buf(struct geo_buf *geo);

/*
 * block until the geo_buf with the given frame_stamp has rendered.
 *
 * This can only be called from outside of the gfx_thread.
 * It should only be called indirectly via gfx_thread_wait_for_geo_buf_stamp.
 */
void rend_wait_for_frame_stamp(unsigned stamp);

/*
 * this implements polygon autosorting on a per-group basis (ideally it would
 * be per-polygon and it would account for polygon splitting, but that is not
 * how it is implemented).  order should be a pre-allocated array with a length
 * of n_groups; after this function returns order will contain the indices of
 * each group ordered from back (minimum Z) to front (maximum Z).
 */
void rend_sort_groups(unsigned *order, struct poly_group const *group,
                      unsigned n_groups);

#endif
