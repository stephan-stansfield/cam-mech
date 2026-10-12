"""Collate a folder of images into a single grid image or an animated GIF.

Opens a GUI (defaulting to results/cams) to choose a folder of images, e.g.
"results/cams/cams_2026-09-17_perc=40, max F=60N/cam_plots". Images are sorted
in natural order (cam_plot_0, cam_plot_1, ..., cam_plot_10, ...) and laid out
row by row. The collage is saved next to the image folder, e.g.
".../cams_2026-09-17_perc=40, max F=60N/cam_plots_collage.png". Choosing the
GIF option instead saves the images as frames of ".../cam_plots_animation.gif".

A folder can also be passed on the command line to skip the folder picker:
    python analysis/collate_images.py "results/cams/<run>/cam_plots" --cols 20 --scale 0.5
    python analysis/collate_images.py "results/cams/<run>/cam_plots" --gif --duration 100
"""

import argparse
import math
import re
import threading
from pathlib import Path

from PIL import Image

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIR = REPO_ROOT / "results" / "cams"
IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}


def natural_key(path):
    """Sort key so that 'cam_plot_2' comes before 'cam_plot_10'."""
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r"(\d+)", path.name)]


def find_images(folder):
    """Return the images in folder, or in its cam_plots/ subfolder if it has one.

    The subfolder takes priority because saved collages live in the run folder itself.
    """
    folder = Path(folder)
    if (folder / "cam_plots").is_dir():
        folder = folder / "cam_plots"
    return sorted((p for p in folder.iterdir() if p.suffix.lower() in IMAGE_EXTS), key=natural_key)


def default_cols(n):
    """Columns for a grid about twice as wide as tall (200 images -> 20 x 10)."""
    return max(1, math.ceil(math.sqrt(2 * n)))


def tile_size(images, scale):
    """Size of every tile: the first image's dimensions times `scale`."""
    with Image.open(images[0]) as first:
        return max(1, round(first.width * scale)), max(1, round(first.height * scale))


def load_tile(path, size, background=(255, 255, 255)):
    """Open an image, resize it to `size` and composite any transparency onto `background`."""
    with Image.open(path) as im:
        im = im.convert("RGBA")
        if im.size != size:
            im = im.resize(size, Image.LANCZOS)
        tile = Image.new("RGB", im.size, background)
        tile.paste(im, mask=im.getchannel("A"))
    return tile


def collate(images, cols, scale=1.0, padding=0, background=(255, 255, 255), progress=None):
    """Paste images into a grid of `cols` columns and return the collage.

    Every tile is sized to the first image's dimensions times `scale`.
    Transparent images are composited onto `background`.
    """
    rows = math.ceil(len(images) / cols)
    tile_w, tile_h = tile_size(images, scale)

    collage = Image.new(
        "RGB",
        (cols * tile_w + (cols + 1) * padding, rows * tile_h + (rows + 1) * padding),
        background,
    )

    for i, path in enumerate(images):
        tile = load_tile(path, (tile_w, tile_h), background)
        row, col = divmod(i, cols)
        collage.paste(tile, (padding + col * (tile_w + padding), padding + row * (tile_h + padding)))
        if progress:
            progress(i + 1, len(images))

    return collage


def make_gif(images, out, scale=1.0, duration=200, background=(255, 255, 255), progress=None):
    """Save the images as an endlessly looping GIF with `duration` ms per frame."""
    size = tile_size(images, scale)
    frames = []
    for i, path in enumerate(images):
        frames.append(load_tile(path, size, background).quantize(colors=256, dither=Image.Dither.NONE))
        if progress:
            progress(i + 1, len(images))
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=duration, loop=0, optimize=True)
    return size


def output_path(folder, suffix):
    """Save beside the image folder so the output isn't picked up on a re-run."""
    folder = Path(folder)
    return folder.parent / f"{folder.name}_{suffix}"


def run(folder, mode="collage", cols=None, scale=1.0, padding=0, duration=200, progress=None):
    """Create a collage or GIF from the images in `folder`. Returns (output path, image count, pixel size)."""
    images = find_images(folder)
    if not images:
        raise ValueError(f"No images found in {folder}")

    if mode == "gif":
        out = output_path(images[0].parent, "animation.gif")
        size = make_gif(images, out, scale, duration, progress=progress)
        return out, len(images), size

    cols = cols or default_cols(len(images))
    collage = collate(images, cols, scale, padding, progress=progress)
    out = output_path(images[0].parent, "collage.png")
    collage.save(out, dpi=(300, 300))
    return out, len(images), collage.size


def launch_gui():
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.title("Collate Images")

    folder_var = tk.StringVar()
    mode_var = tk.StringVar(value="collage")
    cols_var = tk.StringVar()
    scale_var = tk.StringVar(value="0.5")
    padding_var = tk.StringVar(value="0")
    duration_var = tk.StringVar(value="200")
    info_var = tk.StringVar(value="Choose a folder of images.")

    def choose_folder():
        start = DEFAULT_DIR if DEFAULT_DIR.is_dir() else REPO_ROOT
        folder = filedialog.askdirectory(initialdir=start, title="Choose image folder")
        if not folder:
            return
        folder_var.set(folder)
        update_info()

    def update_mode(*_):
        gif = mode_var.get() == "gif"
        for entry in (cols_entry, padding_entry):
            entry.config(state="disabled" if gif else "normal")
        duration_entry.config(state="normal" if gif else "disabled")
        create_btn.config(text="Create GIF" if gif else "Create collage")
        update_info()

    def update_info(*_):
        folder = folder_var.get()
        if not folder:
            return
        images = find_images(folder)
        if not images:
            info_var.set("No images found in this folder.")
            return
        if not cols_var.get():
            cols_var.set(str(default_cols(len(images))))
        try:
            scale = float(scale_var.get())
            with Image.open(images[0]) as im:
                w, h = im.size
            if mode_var.get() == "gif":
                seconds = len(images) * int(duration_var.get()) / 1000
                info_var.set(
                    f"{len(images)} frames ({w}x{h}) -> GIF about "
                    f"{round(w * scale)} x {round(h * scale)} px, {seconds:g} s per loop"
                )
            else:
                cols = int(cols_var.get())
                rows = math.ceil(len(images) / cols)
                info_var.set(
                    f"{len(images)} images ({w}x{h}) -> {cols} x {rows} grid, "
                    f"output about {round(cols * w * scale)} x {round(rows * h * scale)} px"
                )
        except (ValueError, ZeroDivisionError):
            info_var.set(f"{len(images)} images found.")

    def create():
        mode = mode_var.get()
        try:
            cols = int(cols_var.get())
            scale = float(scale_var.get())
            padding = int(padding_var.get())
            duration = int(duration_var.get())
            if cols < 1 or scale <= 0 or padding < 0 or duration < 1:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid input", "Columns, padding and frame duration must be integers and scale must be > 0."
            )
            return

        create_btn.config(state="disabled")

        def progress(done, total):
            root.after(0, lambda: (progress_bar.config(maximum=total, value=done),
                                   info_var.set(f"Processing image {done} of {total}...")))

        def work():
            try:
                out, n, size = run(folder_var.get(), mode, cols, scale, padding, duration, progress)
                msg = f"Saved {n} images ({size[0]}x{size[1]} px) to:\n{out}"
                root.after(0, lambda: (info_var.set(f"Saved {out.name}"), messagebox.showinfo("Done", msg)))
            except Exception as e:
                root.after(0, lambda e=e: messagebox.showerror("Error", str(e)))
            finally:
                root.after(0, lambda: create_btn.config(state="normal"))

        info_var.set("Saving...")
        threading.Thread(target=work, daemon=True).start()

    frame = ttk.Frame(root, padding=12)
    frame.grid(sticky="nsew")
    root.columnconfigure(0, weight=1)
    frame.columnconfigure(1, weight=1)

    ttk.Label(frame, text="Folder:").grid(row=0, column=0, sticky="w")
    ttk.Entry(frame, textvariable=folder_var, width=70).grid(row=0, column=1, sticky="ew", padx=6)
    ttk.Button(frame, text="Browse...", command=choose_folder).grid(row=0, column=2)

    ttk.Label(frame, text="Output:").grid(row=1, column=0, sticky="w", pady=2)
    mode_frame = ttk.Frame(frame)
    mode_frame.grid(row=1, column=1, sticky="w", padx=6)
    ttk.Radiobutton(mode_frame, text="Collage (.png)", variable=mode_var, value="collage").pack(side="left")
    ttk.Radiobutton(mode_frame, text="Animation (.gif)", variable=mode_var, value="gif").pack(side="left", padx=(12, 0))
    mode_var.trace_add("write", update_mode)

    entries = []
    for r, (label, var) in enumerate(
        [
            ("Columns:", cols_var),
            ("Scale (1 = full res):", scale_var),
            ("Padding (px):", padding_var),
            ("Frame duration (ms):", duration_var),
        ],
        start=2,
    ):
        ttk.Label(frame, text=label).grid(row=r, column=0, sticky="w", pady=2)
        entry = ttk.Entry(frame, textvariable=var, width=10)
        entry.grid(row=r, column=1, sticky="w", padx=6)
        entries.append(entry)
        var.trace_add("write", update_info)
    cols_entry, _, padding_entry, duration_entry = entries

    ttk.Label(frame, textvariable=info_var).grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 4))
    progress_bar = ttk.Progressbar(frame, mode="determinate")
    progress_bar.grid(row=7, column=0, columnspan=3, sticky="ew")
    create_btn = ttk.Button(frame, text="Create collage", command=create)
    create_btn.grid(row=8, column=2, pady=(8, 0), sticky="e")

    update_mode()
    root.after(100, choose_folder)
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", nargs="?", help="Folder of images (omit to open the GUI)")
    parser.add_argument("--gif", action="store_true", help="Make an animated GIF instead of a collage")
    parser.add_argument("--cols", type=int, help="Number of collage columns (default: grid about 2:1 wide)")
    parser.add_argument("--scale", type=float, default=0.5, help="Scale applied to each image (default 0.5)")
    parser.add_argument("--padding", type=int, default=0, help="Padding between collage images in pixels")
    parser.add_argument("--duration", type=int, default=200, help="GIF frame duration in ms (default 200)")
    args = parser.parse_args()

    if args.folder is None:
        launch_gui()
        return

    mode = "gif" if args.gif else "collage"
    out, n, size = run(args.folder, mode, args.cols, args.scale, args.padding, args.duration)
    print(f"Saved {n} images ({size[0]}x{size[1]} px) to {out}")


if __name__ == "__main__":
    main()
