import nafuma.auxillary as aux
import nafuma.plotting as btp
import nafuma.eds.io as io

import numpy as np

def show_image(data, options={}):


    default_options = {
        'hide_x_labels': True,
        'hide_y_labels': True,
        'hide_x_ticklabels': True,
        'hide_y_ticklabels': True,
        'hide_x_ticks': True,
        'hide_y_ticks': True,
        'colours': None,
        'brightness': None,
        'show_image': True,
        'resize': None,
        'crop': None,
        'ax': None,
        'fig': None,
    }
    
    options = aux.update_options(options=options, required_options=default_options.keys(), default_options=default_options)

    

    if not isinstance(data['path'], list):
        data['path'] = [data['path']]


    if not 'image' in data.keys():

        data['image'] = [None for _ in range(len(data['path']))]

        if not 'weights' in data.keys():
            data['weights'] = [1.0 for _ in range(len(data['path']))]

        if not options['colours']:
            options['colours'] = [None for _ in range(len(data['path']))]
    
        for i, (path, weight, colour) in enumerate(zip(data['path'], data['weights'], options['colours'])):
            data['image'][i] = io.read_image(path=path, weight=weight, colour=colour, resize=options['resize'], crop=options['crop'])

    
    images = []
    for i, image in enumerate(data['image']):
        images.append(image)
#
    final_image = np.mean(images, axis=0) / 255
    if options['brightness']:
        final_image = io.increase_brightness(final_image, brightness=options['brightness'])

    if len(data['path']) > 1:
        data['image'].append(final_image)


    if options['show_image']:
        if not options['fig'] and not options['ax']:
            fig, ax = btp.prepare_plot(options)
        else:
            fig, ax = options['fig'], options['ax']

        ax.imshow(final_image)
        btp.adjust_plot(fig=fig, ax=ax, options=options)

        return data['image'], fig, ax
    
    else:
        return data['image'], None, None

