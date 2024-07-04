#ifndef TEXTURE_H
#define TEXTURE_H

class Texture{
    public:
        Texture(){}
        Texture(unsigned char *texture,const int texture_width,const int texture_height,const int channels)
        : texture(texture), texture_width(texture_width), texture_height(texture_height), channels(channels){}
        unsigned char* getTexture() const { return texture; }
        const int getTexture_width() const { return texture_width; }
        const int getTexture_height() const { return texture_height; }
        const int getChannels() const { return channels; }
    private:
        unsigned char *texture;
        int texture_width;
        int texture_height;
        int channels;


};

#endif