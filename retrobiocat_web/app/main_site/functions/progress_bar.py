

def set_progress_bar(job, progress, text, classname='progress-bar'):
    job.meta['progress'] = {'progressbar_style': f'width: {progress}%',
                            'progressbar_text': text,
                            'progressbar_classname': classname
                            }
    job.save_meta()

