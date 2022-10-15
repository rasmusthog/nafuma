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



def plot_spectrum(data: dict, options={}):

    default_options = {
        'deconvolutions': None,
        'lines': None,
        'colours': None,
        'xlabel': 'Energy', 'xunit': 'keV', 'xlim': None,
        'ylabel': 'Counts', 'yunit': 'arb. u.', 'ylim': None, 'hide_y_ticklabels': True, 'hide_y_ticks': True,
    }

    options = aux.update_options(options=options, default_options=default_options)

    fig, ax = btp.prepare_plot(options=options)


    spectrum = io.read_spectrum(data['path'])

    if options['deconvolutions']:
        
        deconvolutions = []
        if not isinstance(options['deconvolutions'], list):
            options['deconvolutions'] = [options['deconvolutions']]

        if options['colours'] and (len(options['colours']) != len(options['deconvolutions'])):
            options['colours'] = None

        for deconv in options['deconvolutions']:
            df = io.read_spectrum(deconv)
            deconvolutions.append(df)


    
    spectrum.plot(x='Energy', y='Counts', ax=ax, color='black')

    if options['deconvolutions']:
        if options['colours']:
            for deconv, colour in zip(deconvolutions, options['colours']):
                ax.fill_between(x=deconv['Energy'], y1=deconv['Counts'], y2=0, color=colour, alpha=0.4)
        else:
            for deconv in deconvolutions:
                ax.fill_between(x=deconv['Energy'], y1=deconv['Counts'], y2=0, alpha=0.4)


    if not options['xlim']:
        options['xlim'] = [spectrum['Energy'].min(), spectrum['Energy'].max()]

    if not options['ylim']:
        options['ylim'] = [0, 1.1*spectrum['Counts'].max()]

    if options['lines']:
        for i, (line, energy) in enumerate(options['lines'].items()):
            ax.axvline(x=energy, ls='--', lw=0.5, c='black')
            ax.text(s=line, x=energy, y=(0.9-0.1*i)*options['ylim'][1], fontsize=8)


    
    fig, ax = btp.adjust_plot(fig=fig, ax=ax, options=options)


    return spectrum, fig, ax