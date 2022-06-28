import moviepy.editor as mpy
import numpy as np
files = []
frameList = np.tile(np.arange(0, 150, 1), 3)
for frame in frameList:
    frame = int(frame)
    png = "recordings/frame{}.png".format(frame)
    # matplotlibStyle()
    # fig, ax = plt.subplots(1, 1)
    # fig.set_size_inches(4, 3)
    # overlayColorMap(ax, png, [0, 1], "$\phi$", orientation="vertical")
    # figureName = "colormapped/frame{}.png".format(frame)
    # fig.savefig(figureName, transparent=True, dpi=1500)
    files.append(png)
clip = mpy.ImageSequenceClip(files, fps = 20)
clip.write_gif("test.gif")