import java.nio.file.Path


class FastaUtils {

    static String coreFromHeader(String header) {
        def parts = header.split(/\|/, 2)
        def left = parts[0]
        def right = parts.size() > 1 ? parts[1] : ''
        def slashCount = { String value -> value.findAll('/')?.size() ?: 0 }
        def best = (slashCount(left) >= slashCount(right) ? left : right).trim()
        if (best ==~ /^\d{2}-[A-Za-z0-9]+$/ && parts.size() > 1) {
            return right.trim()
        }
        best
    }

    static String uidFromCore(String core) {
        def normalized = core?.toUpperCase()?.replaceAll(/[^A-Z0-9]/, '') ?: ''
        def digest = java.security.MessageDigest
            .getInstance('SHA-1')
            .digest(normalized.bytes)
            .encodeHex()
            .toString()
            .substring(0, 8)
            .toUpperCase()
        (normalized ?: 'S') + digest
    }

    static List flatten(List values) {
        def flattened = []
        values.each { value ->
            flattened.addAll(value instanceof List || value instanceof Object[] ? value : [value])
        }
        flattened
    }

    static String filename(Object value) {
        if (value instanceof Path) {
            return value.getFileName().toString()
        }
        if (value instanceof File) {
            return value.getName()
        }
        value?.toString()
    }
}
