#!/usr/bin/env python3

import os
import argparse

def TransDecoderTranslateDenovo(file_transcript, 
                                TransDecoderPath = '', 
                                outfolder = '.'
                                ):
    '''
    given a file_transcript, use TransDecoder to translate to proteins and output the result in outfolder
    default, TransDecoder is in $PATH
    outfolder is current directory
    '''
    
    cmd = f'''cd {outfolder} &&\
            
            '''
    
    
    
    
description = '''
given a file_transcript, use TransDecoder to translate to proteins and output the result in outfolder
default, TransDecoder is in $PATH
outfolder is current directory

tested with TransDecoder v5.5.0
'''



